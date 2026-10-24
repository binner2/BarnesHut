/**
 * Benchmark Framework Implementation
 *
 * CPU octree (serial / OpenMP) vs GPU LBVH Barnes-Hut (eager / CUDA Graph)
 * vs GPU direct summation, on several particle distributions.
 *
 * Protocol: warm-up steps are not timed; each measured step is timed
 * separately (CUDA events on the GPU, steady_clock on the CPU) and reported
 * as median and interquartile range. Accuracy is measured on the forces of
 * the common initial state against a double-precision direct sum.
 */

#include "benchmark_framework.cuh"
#include "bh_gpu.cuh"
#include "../tree.h"
#include "../particle.h"
#include "../stdinc.h"

#include <cuda_runtime.h>
#include <chrono>
#include <random>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace barnes_hut {
namespace benchmark {

// ============================================================================
// Utilities
// ============================================================================

const char* distribution_name(Distribution dist) {
    switch (dist) {
        case Distribution::Uniform:    return "Uniform";
        case Distribution::Clustered:  return "Clustered";
        case Distribution::Plummer:    return "Plummer";
        case Distribution::Disk:       return "Disk";
        case Distribution::Shell:      return "Shell";
        case Distribution::DualGalaxy: return "DualGalaxy";
        default: return "Unknown";
    }
}

const char* distribution_regularity(Distribution dist) {
    switch (dist) {
        case Distribution::Uniform:    return "Duzenli (Regular)";
        case Distribution::Disk:       return "Yari-Duzenli (Semi-Regular)";
        case Distribution::Shell:      return "Yari-Duzenli (Semi-Regular)";
        case Distribution::Clustered:  return "Duzensiz (Irregular)";
        case Distribution::Plummer:    return "Duzensiz (Irregular)";
        case Distribution::DualGalaxy: return "Duzensiz (Irregular)";
        default: return "Unknown";
    }
}

static double now_ms() {
    auto t = std::chrono::steady_clock::now();
    return std::chrono::duration<double, std::milli>(t.time_since_epoch()).count();
}

// Linear-interpolated quantile of a sorted sample
static double quantile(const std::vector<double>& sorted, double q) {
    if (sorted.empty()) return 0.0;
    const double pos = q * static_cast<double>(sorted.size() - 1);
    const size_t lo = static_cast<size_t>(pos);
    const size_t hi = std::min(lo + 1, sorted.size() - 1);
    return sorted[lo] + (pos - static_cast<double>(lo)) * (sorted[hi] - sorted[lo]);
}

static TimingResult compute_timing_stats(const std::vector<double>& samples) {
    TimingResult result{};
    if (samples.empty()) return result;

    result.num_samples = static_cast<int>(samples.size());

    double sum = 0.0;
    for (auto s : samples) sum += s;
    result.mean_ms = sum / samples.size();

    double sq_sum = 0.0;
    for (auto s : samples) sq_sum += (s - result.mean_ms) * (s - result.mean_ms);
    result.std_dev_ms = samples.size() > 1 ? std::sqrt(sq_sum / (samples.size() - 1)) : 0.0;

    auto sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    result.min_ms = sorted.front();
    result.max_ms = sorted.back();
    result.median_ms = quantile(sorted, 0.50);
    result.p25_ms = quantile(sorted, 0.25);
    result.p75_ms = quantile(sorted, 0.75);

    return result;
}

// ============================================================================
// Data Generators
// ============================================================================

void generate_distribution(
    Distribution dist, size_t count,
    float* mass, float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    unsigned int seed)
{
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> uni01(0.0f, 1.0f);
    std::normal_distribution<float> normal(0.0f, 1.0f);
    std::uniform_real_distribution<float> mass_dist(1000.0f, 10000.0f);

    for (size_t i = 0; i < count; ++i) {
        mass[i] = mass_dist(rng);
        vel_x[i] = normal(rng) * 0.1f;
        vel_y[i] = normal(rng) * 0.1f;
        vel_z[i] = normal(rng) * 0.1f;
    }

    switch (dist) {
        case Distribution::Uniform: {
            // Küp içinde eşit dağılım [0, 100]^3
            std::uniform_real_distribution<float> pos_dist(0.0f, 100.0f);
            for (size_t i = 0; i < count; ++i) {
                pos_x[i] = pos_dist(rng);
                pos_y[i] = pos_dist(rng);
                pos_z[i] = pos_dist(rng);
            }
            break;
        }

        case Distribution::Clustered: {
            // 5 küme, her kümenin merkezi rastgele, sigma küçük
            const int num_clusters = 5;
            float cx[5], cy[5], cz[5];
            std::uniform_real_distribution<float> center_dist(10.0f, 90.0f);
            for (int c = 0; c < num_clusters; ++c) {
                cx[c] = center_dist(rng);
                cy[c] = center_dist(rng);
                cz[c] = center_dist(rng);
            }

            std::normal_distribution<float> cluster_spread(0.0f, 5.0f);
            for (size_t i = 0; i < count; ++i) {
                int c = static_cast<int>(i % num_clusters);
                pos_x[i] = cx[c] + cluster_spread(rng);
                pos_y[i] = cy[c] + cluster_spread(rng);
                pos_z[i] = cz[c] + cluster_spread(rng);
            }
            break;
        }

        case Distribution::Plummer: {
            // Plummer sphere model (standard astrophysics test)
            // r = a / sqrt(U^(-2/3) - 1) where U ~ Uniform(0,1)
            const float a = 10.0f;  // Plummer radius
            for (size_t i = 0; i < count; ++i) {
                float u = uni01(rng);
                while (u < 1e-6f) u = uni01(rng);  // Avoid singularity
                float r = a / std::sqrt(std::pow(u, -2.0f/3.0f) - 1.0f);

                // Random direction on sphere
                float theta = std::acos(2.0f * uni01(rng) - 1.0f);
                float phi = 2.0f * 3.14159265f * uni01(rng);

                pos_x[i] = 50.0f + r * std::sin(theta) * std::cos(phi);
                pos_y[i] = 50.0f + r * std::sin(theta) * std::sin(phi);
                pos_z[i] = 50.0f + r * std::cos(theta);

                // Plummer velocity (simplified)
                float v_esc = std::sqrt(2.0f) * std::pow(1.0f + r*r/(a*a), -0.25f);
                float v = v_esc * 0.3f;
                vel_x[i] = v * normal(rng);
                vel_y[i] = v * normal(rng);
                vel_z[i] = v * normal(rng);
            }
            break;
        }

        case Distribution::Disk: {
            // Disk galaxy model (exponential radial, thin z)
            const float scale_r = 20.0f;
            const float scale_z = 1.0f;
            for (size_t i = 0; i < count; ++i) {
                // Exponential radial distribution
                float r = -scale_r * std::log(1.0f - uni01(rng) * 0.99f);
                float phi = 2.0f * 3.14159265f * uni01(rng);

                pos_x[i] = 50.0f + r * std::cos(phi);
                pos_y[i] = 50.0f + r * std::sin(phi);
                pos_z[i] = 50.0f + normal(rng) * scale_z;

                // Circular velocity + small random
                float v_circ = std::sqrt(r) * 0.5f;
                vel_x[i] = -v_circ * std::sin(phi) + normal(rng) * 0.05f;
                vel_y[i] =  v_circ * std::cos(phi) + normal(rng) * 0.05f;
                vel_z[i] = normal(rng) * 0.01f;
            }
            break;
        }

        case Distribution::Shell: {
            // Spherical shell
            const float radius = 30.0f;
            const float thickness = 2.0f;
            for (size_t i = 0; i < count; ++i) {
                float r = radius + normal(rng) * thickness;
                float theta = std::acos(2.0f * uni01(rng) - 1.0f);
                float phi = 2.0f * 3.14159265f * uni01(rng);

                pos_x[i] = 50.0f + r * std::sin(theta) * std::cos(phi);
                pos_y[i] = 50.0f + r * std::sin(theta) * std::sin(phi);
                pos_z[i] = 50.0f + r * std::cos(theta);
            }
            break;
        }

        case Distribution::DualGalaxy: {
            // Two colliding Plummer spheres
            const float a = 8.0f;
            const float sep = 40.0f;
            for (size_t i = 0; i < count; ++i) {
                float u = uni01(rng);
                while (u < 1e-6f) u = uni01(rng);
                float r = a / std::sqrt(std::pow(u, -2.0f/3.0f) - 1.0f);

                float theta = std::acos(2.0f * uni01(rng) - 1.0f);
                float phi = 2.0f * 3.14159265f * uni01(rng);

                float offset_x = (i < count / 2) ? -sep/2 : sep/2;
                float drift_vx = (i < count / 2) ? 0.5f : -0.5f;

                pos_x[i] = 50.0f + offset_x + r * std::sin(theta) * std::cos(phi);
                pos_y[i] = 50.0f + r * std::sin(theta) * std::sin(phi);
                pos_z[i] = 50.0f + r * std::cos(theta);

                vel_x[i] = drift_vx + normal(rng) * 0.1f;
                vel_y[i] = normal(rng) * 0.1f;
                vel_z[i] = normal(rng) * 0.1f;
            }
            break;
        }
    }
}

// ============================================================================
// CPU Benchmarks
// ============================================================================

namespace {

#ifdef _OPENMP
// Sets the OpenMP thread count for the lifetime of the object
class ThreadCountScope {
public:
    explicit ThreadCountScope(int threads) : previous_(omp_get_max_threads()) { omp_set_num_threads(threads); }
    ~ThreadCountScope() { omp_set_num_threads(previous_); }
    ThreadCountScope(const ThreadCountScope&) = delete;
    ThreadCountScope& operator=(const ThreadCountScope&) = delete;
private:
    int previous_;
};
#endif

std::vector<Particle> to_particles(const float* mass, const float* px, const float* py, const float* pz,
                                   const float* vx, const float* vy, const float* vz, size_t count) {
    std::vector<Particle> particles;
    particles.reserve(count);
    for (size_t i = 0; i < count; ++i) {
        particles.emplace_back(static_cast<Real>(mass[i]), Vector3D(px[i], py[i], pz[i]),
                               Vector3D(vx[i], vy[i], vz[i]));
    }
    return particles;
}

StepBreakdown benchmark_cpu(
    int threads,
    float* mass, float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count, const BenchmarkConfig& config, int num_steps, int warmup_steps)
{
#ifdef _OPENMP
    ThreadCountScope scope(threads);
#else
    (void)threads;
#endif
    auto particles = to_particles(mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, count);
    BarnesHutTree tree(particles, config.dt, config.theta, config.max_particles_per_leaf);

    for (int s = 0; s < warmup_steps; ++s) {
        tree.simulation_step();
    }

    std::vector<double> total, build, force, integrate;
    for (int s = 0; s < num_steps; ++s) {
        const double t0 = now_ms();
        tree.simulation_step();
        const double t1 = now_ms();
        const auto& st = tree.get_statistics();
        total.push_back(t1 - t0);
        build.push_back(1e3 * (st.time_load + st.time_upward));
        force.push_back(1e3 * st.time_force);
        integrate.push_back(std::max(0.0, (t1 - t0) - 1e3 * (st.time_load + st.time_upward + st.time_force)));
    }

    // Final state back to the caller (positions, velocities and forces)
    for (size_t i = 0; i < count; ++i) {
        force_x[i] = static_cast<float>(particles[i].force()[0]);
        force_y[i] = static_cast<float>(particles[i].force()[1]);
        force_z[i] = static_cast<float>(particles[i].force()[2]);
        pos_x[i] = static_cast<float>(particles[i].position()[0]);
        pos_y[i] = static_cast<float>(particles[i].position()[1]);
        pos_z[i] = static_cast<float>(particles[i].position()[2]);
        vel_x[i] = static_cast<float>(particles[i].velocity()[0]);
        vel_y[i] = static_cast<float>(particles[i].velocity()[1]);
        vel_z[i] = static_cast<float>(particles[i].velocity()[2]);
    }

    StepBreakdown result;
    result.total = compute_timing_stats(total);
    result.tree_build = compute_timing_stats(build);
    result.force_calc = compute_timing_stats(force);
    result.integration = compute_timing_stats(integrate);
    return result;
}

gpu::GpuConfig to_gpu_config(const BenchmarkConfig& config, bool use_graph) {
    gpu::GpuConfig g;
    g.theta = config.theta;
    g.epsilon_squared = config.epsilon_squared;
    g.gravity = config.gravity;
    g.dt = config.dt;
    g.use_graph = use_graph;
    return g;
}

} // namespace

StepBreakdown benchmark_cpu_serial(
    float* mass, float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count, const BenchmarkConfig& config, int num_steps, int warmup_steps)
{
    return benchmark_cpu(1, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                         force_x, force_y, force_z, count, config, num_steps, warmup_steps);
}

StepBreakdown benchmark_cpu_openmp(
    float* mass, float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count, const BenchmarkConfig& config, int num_steps, int warmup_steps)
{
#ifdef _OPENMP
    const int threads = omp_get_max_threads();
#else
    const int threads = 1;
#endif
    return benchmark_cpu(threads, mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z,
                         force_x, force_y, force_z, count, config, num_steps, warmup_steps);
}

// ============================================================================
// GPU Benchmarks
// ============================================================================

bool gpu_available() {
    return gpu::BarnesHutGpu::device_available();
}

void cpu_bh_forces(const float* mass, const float* pos_x, const float* pos_y, const float* pos_z,
                   size_t count, const BenchmarkConfig& config,
                   float* force_x, float* force_y, float* force_z)
{
    const std::vector<float> zero(count, 0.0f);
    auto particles = to_particles(mass, pos_x, pos_y, pos_z, zero.data(), zero.data(), zero.data(), count);
    BarnesHutTree tree(particles, config.dt, config.theta, config.max_particles_per_leaf);
    tree.compute_forces();
    for (size_t i = 0; i < count; ++i) {
        force_x[i] = static_cast<float>(particles[i].force()[0]);
        force_y[i] = static_cast<float>(particles[i].force()[1]);
        force_z[i] = static_cast<float>(particles[i].force()[2]);
    }
}

bool gpu_bh_forces(const float* mass, const float* pos_x, const float* pos_y, const float* pos_z,
                   size_t count, const BenchmarkConfig& config,
                   float* force_x, float* force_y, float* force_z)
{
    if (!gpu_available()) {
        return false;
    }
    const std::vector<float> zero(count, 0.0f);
    gpu::BarnesHutGpu sim(mass, pos_x, pos_y, pos_z, zero.data(), zero.data(), zero.data(), count,
                          to_gpu_config(config, false));
    sim.download(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, force_x, force_y, force_z);
    return true;
}

StepBreakdown benchmark_gpu_direct(
    float* h_mass, float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count, const BenchmarkConfig& config, int num_steps, int warmup_steps)
{
    // Times O(N^2) force evaluations; the state is not advanced, so the
    // returned forces belong to the input positions.
    gpu::BarnesHutGpu sim(h_mass, h_pos_x, h_pos_y, h_pos_z, h_vel_x, h_vel_y, h_vel_z,
                          count, to_gpu_config(config, false));
    for (int s = 0; s < warmup_steps; ++s) {
        sim.compute_forces_direct();
        sim.synchronize();
    }
    std::vector<double> total;
    for (int s = 0; s < num_steps; ++s) {
        sim.compute_forces_direct();
        total.push_back(sim.synchronize().force_ms);
    }
    sim.download(h_pos_x, h_pos_y, h_pos_z, h_vel_x, h_vel_y, h_vel_z, h_force_x, h_force_y, h_force_z);

    StepBreakdown result;
    result.total = compute_timing_stats(total);
    result.force_calc = result.total;
    return result;
}

StepBreakdown benchmark_gpu_bh(
    float* h_mass, float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count, const BenchmarkConfig& config, int num_steps, int warmup_steps, bool use_graph)
{
    gpu::BarnesHutGpu sim(h_mass, h_pos_x, h_pos_y, h_pos_z, h_vel_x, h_vel_y, h_vel_z,
                          count, to_gpu_config(config, use_graph));
    for (int s = 0; s < warmup_steps; ++s) {
        sim.step();
        sim.synchronize();
    }
    std::vector<double> total, build, force, integrate;
    for (int s = 0; s < num_steps; ++s) {
        sim.step();
        const auto t = sim.synchronize();
        total.push_back(t.total_ms);
        if (!use_graph) {
            build.push_back(t.tree_ms);
            force.push_back(t.force_ms);
            integrate.push_back(t.integrate_ms);
        }
    }
    sim.download(h_pos_x, h_pos_y, h_pos_z, h_vel_x, h_vel_y, h_vel_z, h_force_x, h_force_y, h_force_z);
    if (sim.force_stats().stack_overflows > 0) {
        std::cerr << "    warning: GPU traversal stack overflow (forces approximated)\n";
    }

    StepBreakdown result;
    result.total = compute_timing_stats(total);
    result.tree_build = compute_timing_stats(build);
    result.force_calc = compute_timing_stats(force);
    result.integration = compute_timing_stats(integrate);
    return result;
}

// ============================================================================
// Correctness Validation
// ============================================================================

void reference_forces_direct(
    const float* mass, const float* pos_x, const float* pos_y, const float* pos_z,
    size_t count, float epsilon_squared, float gravity,
    const std::vector<size_t>& indices,
    float* force_x, float* force_y, float* force_z)
{
    const long long m = static_cast<long long>(indices.size());
    #ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic, 16)
    #endif
    for (long long k = 0; k < m; ++k) {
        const size_t i = indices[static_cast<size_t>(k)];
        double fx = 0.0, fy = 0.0, fz = 0.0;
        for (size_t j = 0; j < count; ++j) {
            if (j == i) continue;
            const double dx = double(pos_x[j]) - pos_x[i];
            const double dy = double(pos_y[j]) - pos_y[i];
            const double dz = double(pos_z[j]) - pos_z[i];
            const double r2 = dx * dx + dy * dy + dz * dz + epsilon_squared;
            const double s = gravity * double(mass[i]) * mass[j] / (r2 * std::sqrt(r2));
            fx += s * dx;
            fy += s * dy;
            fz += s * dz;
        }
        force_x[k] = static_cast<float>(fx);
        force_y[k] = static_cast<float>(fy);
        force_z[k] = static_cast<float>(fz);
    }
}

CorrectnessResult validate_forces(
    const float* ref_fx, const float* ref_fy, const float* ref_fz,
    const float* test_fx, const float* test_fy, const float* test_fz,
    size_t count, float tolerance)
{
    CorrectnessResult result{};
    double sum_ref_sq = 0.0, sum_diff_sq = 0.0, sum_rel = 0.0;
    size_t valid = 0;
    bool finite = true;

    for (size_t i = 0; i < count; ++i) {
        const double r[3] = {ref_fx[i], ref_fy[i], ref_fz[i]};
        const double t[3] = {test_fx[i], test_fy[i], test_fz[i]};
        if (!std::isfinite(r[0]) || !std::isfinite(r[1]) || !std::isfinite(r[2]) ||
            !std::isfinite(t[0]) || !std::isfinite(t[1]) || !std::isfinite(t[2])) {
            finite = false;
            continue;
        }
        const double ref2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        const double d2 = (r[0] - t[0]) * (r[0] - t[0]) + (r[1] - t[1]) * (r[1] - t[1]) + (r[2] - t[2]) * (r[2] - t[2]);
        sum_ref_sq += ref2;
        sum_diff_sq += d2;
        if (ref2 > 0.0) {
            const double rel = std::sqrt(d2 / ref2);
            result.max_relative_error = std::max(result.max_relative_error, rel);
            sum_rel += rel;
            ++valid;
        }
    }

    result.mean_relative_error = valid ? sum_rel / valid : 0.0;
    result.l2_norm_error = sum_ref_sq > 0.0 ? std::sqrt(sum_diff_sq / sum_ref_sq) : std::sqrt(sum_diff_sq);
    if (!finite) {
        result.max_relative_error = std::numeric_limits<double>::infinity();
    }
    // Global relative L2 error ||F - F_ref|| / ||F_ref||; per-particle maxima
    // are dominated by particles with near-zero net force and are reported
    // but not used as the pass criterion.
    result.passed = finite && result.l2_norm_error <= tolerance;
    return result;
}

// ============================================================================
// Full Benchmark Suite
// ============================================================================

namespace {

// Particles used for force validation: all of them up to the configured
// limit, otherwise an evenly spaced sample (direct reference is O(N * S)).
std::vector<size_t> validation_indices(size_t n, size_t limit) {
    std::vector<size_t> idx;
    if (n <= limit) {
        idx.resize(n);
        std::iota(idx.begin(), idx.end(), size_t{0});
    } else {
        for (size_t k = 0; k < limit; ++k) idx.push_back(k * n / limit);
    }
    return idx;
}

CorrectnessResult compare_subset(const std::vector<size_t>& idx,
                                 const std::vector<float>& rfx, const std::vector<float>& rfy,
                                 const std::vector<float>& rfz,
                                 const std::vector<float>& fx, const std::vector<float>& fy,
                                 const std::vector<float>& fz, float tolerance) {
    std::vector<float> tx(idx.size()), ty(idx.size()), tz(idx.size());
    for (size_t k = 0; k < idx.size(); ++k) {
        tx[k] = fx[idx[k]]; ty[k] = fy[idx[k]]; tz[k] = fz[idx[k]];
    }
    return validate_forces(rfx.data(), rfy.data(), rfz.data(), tx.data(), ty.data(), tz.data(),
                           idx.size(), tolerance);
}

BenchmarkResult make_result(const std::string& method, Distribution dist, size_t n,
                            const StepBreakdown& timing, const CorrectnessResult& correctness,
                            double serial_median) {
    BenchmarkResult r;
    r.method_name = method;
    r.distribution = dist;
    r.particle_count = n;
    r.timing = timing;
    r.correctness = correctness;
    const double t = std::max(timing.total.median_ms, 1e-6);
    r.speedup_vs_serial = serial_median > 0.0 ? serial_median / t : 1.0;
    r.particles_per_second = static_cast<double>(n) / (t / 1000.0);
    r.gflops = 0.0;
    return r;
}

void print_row(const BenchmarkResult& r) {
    std::cout << "    " << std::left << std::setw(16) << r.method_name << std::right << std::fixed
              << std::setprecision(3) << std::setw(10) << r.timing.total.median_ms << " ms  [IQR "
              << (r.timing.total.p75_ms - r.timing.total.p25_ms) << "]  "
              << std::setprecision(2) << std::setw(8) << r.speedup_vs_serial << "x  L2err "
              << std::scientific << std::setprecision(2) << r.correctness.l2_norm_error
              << (r.correctness.passed ? "" : "  [ACCURACY FAIL]") << std::fixed << "\n";
}

} // namespace

std::vector<BenchmarkResult> run_benchmark_suite(const BenchmarkConfig& config) {
    std::vector<BenchmarkResult> all_results;

    std::cout << "\n"
              << "================================================================\n"
              << "  Barnes-Hut Benchmark Suite (CPU octree vs GPU LBVH)\n"
              << "  warmup=" << config.warmup_steps << " measured=" << config.measured_steps
              << " theta=" << config.theta << " dt=" << config.dt << "\n"
              << "================================================================\n\n";

    const bool has_gpu = gpu_available();
    if (has_gpu) {
        int device = 0;
        cudaGetDevice(&device);
        cudaDeviceProp prop{};
        cudaGetDeviceProperties(&prop, device);
        int runtime = 0, driver = 0;
        cudaRuntimeGetVersion(&runtime);
        cudaDriverGetVersion(&driver);
        std::cout << "GPU: " << prop.name << " (CC " << prop.major << "." << prop.minor << ", "
                  << prop.multiProcessorCount << " SMs, " << (prop.totalGlobalMem >> 20) << " MB)"
                  << " runtime " << runtime << " driver " << driver << "\n";
    } else {
        std::cout << "GPU: not available - GPU benchmarks are skipped\n";
    }
#ifdef _OPENMP
    std::cout << "OpenMP: " << omp_get_max_threads() << " threads\n\n";
#else
    std::cout << "OpenMP: not available\n\n";
#endif

    for (auto dist : config.distributions) {
        std::cout << "--- Distribution: " << distribution_name(dist)
                  << " (" << distribution_regularity(dist) << ") ---\n\n";

        for (auto N : config.particle_counts) {
            std::cout << "  N = " << N << ":\n";

            std::vector<float> mass(N), px0(N), py0(N), pz0(N), vx0(N), vy0(N), vz0(N);
            generate_distribution(dist, N, mass.data(), px0.data(), py0.data(), pz0.data(),
                                  vx0.data(), vy0.data(), vz0.data(), 42);

            // Accuracy is judged on forces at the common initial state
            // against a double-precision direct sum.
            const auto idx = validation_indices(N, config.validation_sample);
            std::vector<float> rfx(idx.size()), rfy(idx.size()), rfz(idx.size());
            reference_forces_direct(mass.data(), px0.data(), py0.data(), pz0.data(), N,
                                    config.epsilon_squared, config.gravity, idx,
                                    rfx.data(), rfy.data(), rfz.data());

            std::vector<float> fx(N), fy(N), fz(N);
            CorrectnessResult cpu_accuracy;
            {
                auto particles = to_particles(mass.data(), px0.data(), py0.data(), pz0.data(),
                                              vx0.data(), vy0.data(), vz0.data(), N);
                BarnesHutTree tree(particles, config.dt, config.theta, config.max_particles_per_leaf);
                tree.compute_forces();
                for (size_t i = 0; i < N; ++i) {
                    fx[i] = float(particles[i].force()[0]);
                    fy[i] = float(particles[i].force()[1]);
                    fz[i] = float(particles[i].force()[2]);
                }
                cpu_accuracy = compare_subset(idx, rfx, rfy, rfz, fx, fy, fz, config.max_l2_error);
            }

            auto run = [&](auto&& bench) {
                auto px = px0, py = py0, pz = pz0, vx = vx0, vy = vy0, vz = vz0;
                return bench(px, py, pz, vx, vy, vz);
            };

            const auto serial = run([&](auto& px, auto& py, auto& pz, auto& vx, auto& vy, auto& vz) {
                return benchmark_cpu_serial(mass.data(), px.data(), py.data(), pz.data(), vx.data(),
                                            vy.data(), vz.data(), fx.data(), fy.data(), fz.data(), N,
                                            config, config.measured_steps, config.warmup_steps);
            });
            all_results.push_back(make_result("CPU_Serial", dist, N, serial, cpu_accuracy, 0.0));
            const double serial_median = serial.total.median_ms;
            print_row(all_results.back());

#ifdef _OPENMP
            const auto omp = run([&](auto& px, auto& py, auto& pz, auto& vx, auto& vy, auto& vz) {
                return benchmark_cpu_openmp(mass.data(), px.data(), py.data(), pz.data(), vx.data(),
                                            vy.data(), vz.data(), fx.data(), fy.data(), fz.data(), N,
                                            config, config.measured_steps, config.warmup_steps);
            });
            all_results.push_back(make_result("CPU_OpenMP", dist, N, omp, cpu_accuracy, serial_median));
            print_row(all_results.back());
#endif

            if (!has_gpu) {
                std::cout << "\n";
                continue;
            }

            try {
                // GPU direct: force evaluations only, state unchanged
                auto px = px0, py = py0, pz = pz0, vx = vx0, vy = vy0, vz = vz0;
                const auto direct = benchmark_gpu_direct(mass.data(), px.data(), py.data(), pz.data(),
                                                         vx.data(), vy.data(), vz.data(), fx.data(),
                                                         fy.data(), fz.data(), N, config,
                                                         config.measured_steps, config.warmup_steps);
                all_results.push_back(make_result("GPU_Direct", dist, N, direct,
                                                  compare_subset(idx, rfx, rfy, rfz, fx, fy, fz, 1e-4f),
                                                  serial_median));
                print_row(all_results.back());

                // GPU Barnes-Hut accuracy: forces computed by the constructor
                CorrectnessResult gpu_accuracy;
                {
                    px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;
                    gpu::BarnesHutGpu sim(mass.data(), px.data(), py.data(), pz.data(), vx.data(),
                                          vy.data(), vz.data(), N, to_gpu_config(config, false));
                    sim.download(nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,
                                 fx.data(), fy.data(), fz.data());
                    gpu_accuracy = compare_subset(idx, rfx, rfy, rfz, fx, fy, fz, config.max_l2_error);
                }

                for (bool graph : {false, true}) {
                    const auto bh = run([&](auto& px_, auto& py_, auto& pz_, auto& vx_, auto& vy_, auto& vz_) {
                        return benchmark_gpu_bh(mass.data(), px_.data(), py_.data(), pz_.data(),
                                                vx_.data(), vy_.data(), vz_.data(), fx.data(),
                                                fy.data(), fz.data(), N, config,
                                                config.measured_steps, config.warmup_steps, graph);
                    });
                    all_results.push_back(make_result(graph ? "GPU_BH_Graph" : "GPU_BH", dist, N, bh,
                                                      gpu_accuracy, serial_median));
                    print_row(all_results.back());
                }
            } catch (const std::exception& e) {
                std::cerr << "    GPU benchmark failed: " << e.what() << "\n";
            }

            std::cout << "\n";
        }
    }

    return all_results;
}

// ============================================================================
// Output Functions
// ============================================================================

void write_results_csv(
    const std::vector<BenchmarkResult>& results,
    const std::string& filename)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Cannot open " << filename << " for writing\n";
        return;
    }

    // Header
    out << "Method,Distribution,Regularity,ParticleCount,"
        << "MedianTime_ms,P25_ms,P75_ms,MeanTime_ms,StdDev_ms,MinTime_ms,MaxTime_ms,Samples,"
        << "SpeedupVsSerial,ParticlesPerSecond,"
        << "MaxRelativeError,MeanRelativeError,L2NormError,CorrectnessPassed,"
        << "TreeBuildMedian_ms,ForceCalcMedian_ms,IntegrationMedian_ms\n";

    for (const auto& r : results) {
        out << r.method_name << ","
            << distribution_name(r.distribution) << ","
            << distribution_regularity(r.distribution) << ","
            << r.particle_count << ","
            << r.timing.total.median_ms << ","
            << r.timing.total.p25_ms << ","
            << r.timing.total.p75_ms << ","
            << r.timing.total.mean_ms << ","
            << r.timing.total.std_dev_ms << ","
            << r.timing.total.min_ms << ","
            << r.timing.total.max_ms << ","
            << r.timing.total.num_samples << ","
            << r.speedup_vs_serial << ","
            << r.particles_per_second << ","
            << r.correctness.max_relative_error << ","
            << r.correctness.mean_relative_error << ","
            << r.correctness.l2_norm_error << ","
            << (r.correctness.passed ? "PASS" : "FAIL") << ","
            << r.timing.tree_build.median_ms << ","
            << r.timing.force_calc.median_ms << ","
            << r.timing.integration.median_ms
            << "\n";
    }

    std::cout << "\nResults written to " << filename << "\n";
}

void print_results_table(const std::vector<BenchmarkResult>& results) {
    std::cout << "\n"
              << "================================================================\n"
              << "  SONUC TABLOSU (Results Table)\n"
              << "================================================================\n\n";

    // Header
    std::cout << std::left
              << std::setw(22) << "Method"
              << std::setw(12) << "Dist"
              << std::setw(10) << "N"
              << std::right
              << std::setw(12) << "Median(ms)"
              << std::setw(10) << "IQR"
              << std::setw(10) << "Speedup"
              << std::setw(15) << "P/s"
              << std::setw(12) << "L2Err"
              << std::setw(8) << "Pass"
              << "\n";
    std::cout << std::string(109, '-') << "\n";

    for (const auto& r : results) {
        std::cout << std::left
                  << std::setw(22) << r.method_name
                  << std::setw(12) << distribution_name(r.distribution)
                  << std::setw(10) << r.particle_count
                  << std::right << std::fixed
                  << std::setw(12) << std::setprecision(2) << r.timing.total.median_ms
                  << std::setw(10) << std::setprecision(2) << (r.timing.total.p75_ms - r.timing.total.p25_ms)
                  << std::setw(10) << std::setprecision(2) << r.speedup_vs_serial
                  << std::setw(15) << std::setprecision(0) << r.particles_per_second
                  << std::setw(12) << std::scientific << std::setprecision(2) << r.correctness.l2_norm_error
                  << std::setw(8) << (r.correctness.passed ? "OK" : "FAIL")
                  << "\n";
    }
}

void print_distribution_analysis(const std::vector<BenchmarkResult>& results) {
    std::cout << "\n"
              << "================================================================\n"
              << "  DAGILIM ANALIZI (Distribution Analysis)\n"
              << "  Duzenli vs Duzensiz Veri Dagilimi Etki Raporu\n"
              << "================================================================\n\n";

    // Group results by method and particle count
    // Compare uniform (regular) vs others (irregular)
    std::cout << "Bu analiz, veri dagılımının (düzenli/düzensiz) farklı\n"
              << "hesaplama yöntemlerindeki performans etkisini gösterir.\n\n";

    // Find all methods
    std::vector<std::string> methods;
    for (const auto& r : results) {
        if (std::find(methods.begin(), methods.end(), r.method_name) == methods.end()) {
            methods.push_back(r.method_name);
        }
    }

    for (const auto& method : methods) {
        std::cout << "--- " << method << " ---\n\n";

        std::cout << std::left
                  << std::setw(12) << "N"
                  << std::setw(14) << "Uniform(ms)"
                  << std::setw(14) << "Clustered"
                  << std::setw(14) << "Plummer"
                  << std::setw(14) << "Disk"
                  << std::setw(14) << "Impact"
                  << "\n";
        std::cout << std::string(82, '-') << "\n";

        for (auto N : {(size_t)1000, (size_t)5000, (size_t)10000, (size_t)50000, (size_t)100000}) {
            float uniform_time = -1, clustered_time = -1, plummer_time = -1, disk_time = -1;

            for (const auto& r : results) {
                if (r.method_name != method || r.particle_count != N) continue;
                switch (r.distribution) {
                    case Distribution::Uniform:   uniform_time = r.timing.total.median_ms; break;
                    case Distribution::Clustered:  clustered_time = r.timing.total.median_ms; break;
                    case Distribution::Plummer:    plummer_time = r.timing.total.median_ms; break;
                    case Distribution::Disk:       disk_time = r.timing.total.median_ms; break;
                    default: break;
                }
            }

            if (uniform_time < 0) continue;

            // Impact = max irregular time / uniform time
            float max_irregular = std::max({clustered_time, plummer_time, disk_time});
            float impact = (uniform_time > 0) ? max_irregular / uniform_time : 0;

            std::cout << std::left << std::setw(12) << N << std::fixed << std::setprecision(2)
                      << std::setw(14) << uniform_time
                      << std::setw(14) << (clustered_time >= 0 ? std::to_string((int)(clustered_time*100)/100.0).substr(0,8) : "N/A")
                      << std::setw(14) << (plummer_time >= 0 ? std::to_string((int)(plummer_time*100)/100.0).substr(0,8) : "N/A")
                      << std::setw(14) << (disk_time >= 0 ? std::to_string((int)(disk_time*100)/100.0).substr(0,8) : "N/A")
                      << std::setw(14) << (std::to_string((int)(impact*100)/100.0).substr(0,5) + "x")
                      << "\n";
        }
        std::cout << "\n";
    }

    std::cout << "YORUM:\n"
              << "- Impact > 1.0: Duzensiz dagılım daha yavas (beklenen)\n"
              << "- CPU icin: Duzensiz dagılım agac derinligini artırır → daha fazla traversal\n"
              << "- GPU icin: Duzensiz dagılım warp divergence yaratır (hipotez; ncu ile dogrulanmali)\n"
              << "\n";
}

} // namespace benchmark
} // namespace barnes_hut
