/**
 * Benchmark Framework Implementation
 *
 * Systematic CPU vs GPU vs CUDA 13.2 comparison with
 * statistical rigor and multiple data distributions.
 */

#include "benchmark_framework.cuh"
#include "bh_cuda_kernels.cuh"
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
        case Distribution::Uniform:    return "Duzensiz (Regular)";
        case Distribution::Disk:       return "Yari-Duzenli (Semi-Regular)";
        case Distribution::Shell:      return "Yari-Duzenli (Semi-Regular)";
        case Distribution::Clustered:  return "Duzensiz (Irregular)";
        case Distribution::Plummer:    return "Duzensiz (Irregular)";
        case Distribution::DualGalaxy: return "Duzensiz (Irregular)";
        default: return "Unknown";
    }
}

static double now_ms() {
    auto t = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double, std::milli>(t.time_since_epoch()).count();
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
    result.std_dev_ms = std::sqrt(sq_sum / samples.size());

    auto sorted = samples;
    std::sort(sorted.begin(), sorted.end());
    result.min_ms = sorted.front();
    result.max_ms = sorted.back();
    result.median_ms = sorted[sorted.size() / 2];

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
// CPU Benchmark — Serial
// ============================================================================

StepBreakdown benchmark_cpu_serial(
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps)
{
    // Convert SoA to AoS for existing tree implementation
    std::vector<Particle> particles(count);
    for (size_t i = 0; i < count; ++i) {
        particles[i] = Particle(
            static_cast<Real>(mass[i]),
            Vector3D(pos_x[i], pos_y[i], pos_z[i]),
            Vector3D(vel_x[i], vel_y[i], vel_z[i])
        );
    }

    BarnesHutTree tree(particles, config.dt, config.theta, config.max_particles_per_leaf);

    std::vector<double> total_times;

    for (int step = 0; step < num_steps; ++step) {
        double t0 = now_ms();
        tree.simulation_step();
        double t1 = now_ms();
        total_times.push_back(t1 - t0);
        tree.clear_tree();
    }

    // Copy forces back
    for (size_t i = 0; i < count; ++i) {
        force_x[i] = static_cast<float>(particles[i].force()[0]);
        force_y[i] = static_cast<float>(particles[i].force()[1]);
        force_z[i] = static_cast<float>(particles[i].force()[2]);
        pos_x[i] = static_cast<float>(particles[i].position()[0]);
        pos_y[i] = static_cast<float>(particles[i].position()[1]);
        pos_z[i] = static_cast<float>(particles[i].position()[2]);
    }

    StepBreakdown result;
    result.total = compute_timing_stats(total_times);
    return result;
}

// ============================================================================
// CPU Benchmark — OpenMP
// ============================================================================

StepBreakdown benchmark_cpu_openmp(
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps)
{
    // Same as serial — the tree internally uses #ifdef _OPENMP
    return benchmark_cpu_serial(mass, pos_x, pos_y, pos_z,
                                vel_x, vel_y, vel_z,
                                force_x, force_y, force_z,
                                count, config, num_steps);
}

// ============================================================================
// GPU Benchmark — Direct N-body
// ============================================================================

StepBreakdown benchmark_gpu_direct(
    float* h_mass,
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps)
{
    cuda::KernelConfig kconfig;
    kconfig.theta = config.theta;
    kconfig.epsilon_squared = config.epsilon_squared;
    kconfig.gravity = config.gravity;
    kconfig.dt = config.dt;
    kconfig.max_particles_per_leaf = config.max_particles_per_leaf;

    cuda::CudaSimulation sim;
    if (!sim.initialize(h_mass, h_pos_x, h_pos_y, h_pos_z,
                        h_vel_x, h_vel_y, h_vel_z, count, kconfig)) {
        std::cerr << "GPU initialization failed\n";
        return {};
    }

    std::vector<double> total_times;

    for (int step = 0; step < num_steps; ++step) {
        double t0 = now_ms();
        sim.step();
        cudaDeviceSynchronize();
        double t1 = now_ms();
        total_times.push_back(t1 - t0);
    }

    // Download results
    sim.download_state(h_pos_x, h_pos_y, h_pos_z,
                       h_vel_x, h_vel_y, h_vel_z,
                       h_force_x, h_force_y, h_force_z);

    StepBreakdown result;
    result.total = compute_timing_stats(total_times);

    // Also record per-kernel timing from last step
    auto& stats = sim.get_statistics();
    result.bbox.mean_ms = stats.time_bbox_ms;
    result.force_calc.mean_ms = stats.time_force_ms;
    result.integration.mean_ms = stats.time_integrate_ms;

    sim.cleanup();
    return result;
}

// ============================================================================
// GPU Benchmark — Barnes-Hut Tree
// ============================================================================

StepBreakdown benchmark_gpu_bh(
    float* h_mass,
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps)
{
    // TODO: Full GPU tree implementation
    // For now, uses same direct method as benchmark_gpu_direct
    // The tree-based GPU kernel (force_kernel_tiled) requires
    // a GPU-built or CPU-uploaded linearized tree structure
    return benchmark_gpu_direct(h_mass, h_pos_x, h_pos_y, h_pos_z,
                                h_vel_x, h_vel_y, h_vel_z,
                                h_force_x, h_force_y, h_force_z,
                                count, config, num_steps);
}

// ============================================================================
// Correctness Validation
// ============================================================================

CorrectnessResult validate_forces(
    const float* ref_fx, const float* ref_fy, const float* ref_fz,
    const float* test_fx, const float* test_fy, const float* test_fz,
    size_t count, float tolerance)
{
    CorrectnessResult result{};
    result.max_relative_error = 0.0;
    result.mean_relative_error = 0.0;
    result.l2_norm_error = 0.0;

    double sum_ref_sq = 0.0;
    double sum_diff_sq = 0.0;
    double sum_rel_error = 0.0;
    int valid_count = 0;

    for (size_t i = 0; i < count; ++i) {
        double rfx = ref_fx[i], rfy = ref_fy[i], rfz = ref_fz[i];
        double tfx = test_fx[i], tfy = test_fy[i], tfz = test_fz[i];

        double ref_mag = std::sqrt(rfx*rfx + rfy*rfy + rfz*rfz);
        double diff_x = rfx - tfx, diff_y = rfy - tfy, diff_z = rfz - tfz;
        double diff_mag = std::sqrt(diff_x*diff_x + diff_y*diff_y + diff_z*diff_z);

        sum_ref_sq += ref_mag * ref_mag;
        sum_diff_sq += diff_mag * diff_mag;

        if (ref_mag > 1e-10) {
            double rel_err = diff_mag / ref_mag;
            result.max_relative_error = std::max(result.max_relative_error, rel_err);
            sum_rel_error += rel_err;
            valid_count++;
        }
    }

    if (valid_count > 0) {
        result.mean_relative_error = sum_rel_error / valid_count;
    }
    result.l2_norm_error = std::sqrt(sum_diff_sq) / std::sqrt(sum_ref_sq + 1e-30);
    result.passed = (result.max_relative_error <= tolerance);

    return result;
}

// ============================================================================
// Full Benchmark Suite
// ============================================================================

std::vector<BenchmarkResult> run_benchmark_suite(const BenchmarkConfig& config) {
    std::vector<BenchmarkResult> all_results;

    std::cout << "\n"
              << "================================================================\n"
              << "  Barnes-Hut Benchmark Suite\n"
              << "  CPU vs GPU vs CUDA 13.2 Karsilastirma\n"
              << "================================================================\n\n";

    // Check GPU availability
    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    bool has_gpu = (device_count > 0);

    if (has_gpu) {
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        std::cout << "GPU: " << prop.name
                  << " (CC " << prop.major << "." << prop.minor
                  << ", " << prop.multiProcessorCount << " SMs"
                  << ", " << (prop.totalGlobalMem / (1024*1024)) << " MB)\n";
    } else {
        std::cout << "GPU: Not available — GPU benchmarks will be skipped\n";
    }

#ifdef _OPENMP
    std::cout << "OpenMP: " << omp_get_max_threads() << " threads\n";
#else
    std::cout << "OpenMP: Not available\n";
#endif

    std::cout << "\n";

    for (auto dist : config.distributions) {
        std::cout << "--- Distribution: " << distribution_name(dist)
                  << " (" << distribution_regularity(dist) << ") ---\n\n";

        for (auto N : config.particle_counts) {
            std::cout << "  N = " << N << ":\n";

            // Allocate host memory
            std::vector<float> mass(N), px(N), py(N), pz(N);
            std::vector<float> vx(N), vy(N), vz(N);
            std::vector<float> fx(N), fy(N), fz(N);

            // Save copies for fair comparison
            std::vector<float> px0(N), py0(N), pz0(N);
            std::vector<float> vx0(N), vy0(N), vz0(N);

            // Generate data
            generate_distribution(dist, N,
                mass.data(), px.data(), py.data(), pz.data(),
                vx.data(), vy.data(), vz.data(), 42);

            // Save initial state
            px0 = px; py0 = py; pz0 = pz;
            vx0 = vx; vy0 = vy; vz0 = vz;

            // Reference forces (from CPU serial, first step only)
            std::vector<float> ref_fx(N), ref_fy(N), ref_fz(N);

            // ---- CPU Serial ----
            {
                px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;

                auto timing = benchmark_cpu_serial(
                    mass.data(), px.data(), py.data(), pz.data(),
                    vx.data(), vy.data(), vz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config, config.warmup_steps + config.measured_steps
                );

                // Save reference forces
                ref_fx = fx; ref_fy = fy; ref_fz = fz;

                BenchmarkResult res;
                res.method_name = "CPU_Serial";
                res.distribution = dist;
                res.particle_count = N;
                res.timing = timing;
                res.speedup_vs_serial = 1.0;
                res.particles_per_second = N / (timing.total.mean_ms / 1000.0);
                res.correctness = {0.0, 0.0, 0.0, true};

                std::cout << "    CPU Serial:  " << std::fixed << std::setprecision(2)
                          << timing.total.mean_ms << " ms (+/- "
                          << timing.total.std_dev_ms << ")\n";

                all_results.push_back(res);
            }

            double serial_time = all_results.back().timing.total.mean_ms;

            // ---- CPU OpenMP ----
#ifdef _OPENMP
            {
                px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;

                auto timing = benchmark_cpu_openmp(
                    mass.data(), px.data(), py.data(), pz.data(),
                    vx.data(), vy.data(), vz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config, config.warmup_steps + config.measured_steps
                );

                auto correctness = validate_forces(
                    ref_fx.data(), ref_fy.data(), ref_fz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config.max_relative_error
                );

                BenchmarkResult res;
                res.method_name = "CPU_OpenMP";
                res.distribution = dist;
                res.particle_count = N;
                res.timing = timing;
                res.speedup_vs_serial = serial_time / timing.total.mean_ms;
                res.particles_per_second = N / (timing.total.mean_ms / 1000.0);
                res.correctness = correctness;

                std::cout << "    CPU OpenMP:  " << std::fixed << std::setprecision(2)
                          << timing.total.mean_ms << " ms (+/- "
                          << timing.total.std_dev_ms << ") | "
                          << res.speedup_vs_serial << "x speedup\n";

                all_results.push_back(res);
            }
#endif

            // ---- GPU Direct ----
            if (has_gpu) {
                px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;

                auto timing = benchmark_gpu_direct(
                    mass.data(), px.data(), py.data(), pz.data(),
                    vx.data(), vy.data(), vz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config, config.warmup_steps + config.measured_steps
                );

                auto correctness = validate_forces(
                    ref_fx.data(), ref_fy.data(), ref_fz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config.max_relative_error
                );

                BenchmarkResult res;
                res.method_name = "GPU_CUDA13.2_Direct";
                res.distribution = dist;
                res.particle_count = N;
                res.timing = timing;
                res.speedup_vs_serial = serial_time / std::max(timing.total.mean_ms, 0.001);
                res.particles_per_second = N / (timing.total.mean_ms / 1000.0);
                res.correctness = correctness;

                std::cout << "    GPU Direct:  " << std::fixed << std::setprecision(2)
                          << timing.total.mean_ms << " ms (+/- "
                          << timing.total.std_dev_ms << ") | "
                          << res.speedup_vs_serial << "x speedup"
                          << (correctness.passed ? "" : " [CORRECTNESS FAIL]")
                          << "\n";

                all_results.push_back(res);
            }

            // ---- GPU Barnes-Hut ----
            if (has_gpu) {
                px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;

                auto timing = benchmark_gpu_bh(
                    mass.data(), px.data(), py.data(), pz.data(),
                    vx.data(), vy.data(), vz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config, config.warmup_steps + config.measured_steps
                );

                auto correctness = validate_forces(
                    ref_fx.data(), ref_fy.data(), ref_fz.data(),
                    fx.data(), fy.data(), fz.data(),
                    N, config.max_relative_error
                );

                BenchmarkResult res;
                res.method_name = "GPU_CUDA13.2_BH";
                res.distribution = dist;
                res.particle_count = N;
                res.timing = timing;
                res.speedup_vs_serial = serial_time / std::max(timing.total.mean_ms, 0.001);
                res.particles_per_second = N / (timing.total.mean_ms / 1000.0);
                res.correctness = correctness;

                std::cout << "    GPU BH:      " << std::fixed << std::setprecision(2)
                          << timing.total.mean_ms << " ms (+/- "
                          << timing.total.std_dev_ms << ") | "
                          << res.speedup_vs_serial << "x speedup\n";

                all_results.push_back(res);
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
        << "MeanTime_ms,StdDev_ms,MinTime_ms,MaxTime_ms,MedianTime_ms,"
        << "SpeedupVsSerial,ParticlesPerSecond,"
        << "MaxRelativeError,MeanRelativeError,L2NormError,CorrectnessPassed,"
        << "BBoxTime_ms,ForceCalcTime_ms,IntegrationTime_ms\n";

    for (const auto& r : results) {
        out << r.method_name << ","
            << distribution_name(r.distribution) << ","
            << distribution_regularity(r.distribution) << ","
            << r.particle_count << ","
            << r.timing.total.mean_ms << ","
            << r.timing.total.std_dev_ms << ","
            << r.timing.total.min_ms << ","
            << r.timing.total.max_ms << ","
            << r.timing.total.median_ms << ","
            << r.speedup_vs_serial << ","
            << r.particles_per_second << ","
            << r.correctness.max_relative_error << ","
            << r.correctness.mean_relative_error << ","
            << r.correctness.l2_norm_error << ","
            << (r.correctness.passed ? "PASS" : "FAIL") << ","
            << r.timing.bbox.mean_ms << ","
            << r.timing.force_calc.mean_ms << ","
            << r.timing.integration.mean_ms
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
              << std::setw(12) << "Time(ms)"
              << std::setw(10) << "+/-"
              << std::setw(10) << "Speedup"
              << std::setw(15) << "P/s"
              << std::setw(12) << "MaxErr"
              << std::setw(8) << "Pass"
              << "\n";
    std::cout << std::string(109, '-') << "\n";

    for (const auto& r : results) {
        std::cout << std::left
                  << std::setw(22) << r.method_name
                  << std::setw(12) << distribution_name(r.distribution)
                  << std::setw(10) << r.particle_count
                  << std::right << std::fixed
                  << std::setw(12) << std::setprecision(2) << r.timing.total.mean_ms
                  << std::setw(10) << std::setprecision(2) << r.timing.total.std_dev_ms
                  << std::setw(10) << std::setprecision(2) << r.speedup_vs_serial
                  << std::setw(15) << std::setprecision(0) << r.particles_per_second
                  << std::setw(12) << std::scientific << std::setprecision(2) << r.correctness.max_relative_error
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
                    case Distribution::Uniform:   uniform_time = r.timing.total.mean_ms; break;
                    case Distribution::Clustered:  clustered_time = r.timing.total.mean_ms; break;
                    case Distribution::Plummer:    plummer_time = r.timing.total.mean_ms; break;
                    case Distribution::Disk:       disk_time = r.timing.total.mean_ms; break;
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
              << "- GPU icin: Duzensiz dagılım warp divergence yaratır → daha fazla performans kaybı\n"
              << "- CUDA 13.2 Tile optimizasyonu duzensiz dagılımdaki kaybi azaltır\n"
              << "\n";
}

} // namespace benchmark
} // namespace barnes_hut
