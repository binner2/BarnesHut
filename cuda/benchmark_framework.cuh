#pragma once

/**
 * Benchmark Framework for Barnes-Hut Simulation
 *
 * Systematic comparison of:
 * 1. CPU Serial (octree, 1 thread)
 * 2. CPU OpenMP (octree, all threads)
 * 3. GPU direct N-body (O(N^2), force evaluation only)
 * 4. GPU Barnes-Hut (LBVH), eager launches and CUDA Graph replay
 *
 * Features:
 * - Untimed warm-up, per-step samples, median + IQR
 * - Data distribution variants: uniform, clustered, Plummer, disk
 * - Scaling analysis: varying particle counts
 * - Accuracy: global relative L2 force error vs double direct sum
 */

#include <string>
#include <vector>
#include <cstddef>
#include <functional>

namespace barnes_hut {
namespace benchmark {

// ============================================================================
// Data Distribution Types
// ============================================================================

enum class Distribution {
    Uniform,          // Düzenli: Küp içinde eşit dağılım
    Clustered,        // Düzensiz: Birden fazla küme etrafında yoğunlaşmış
    Plummer,          // Düzensiz: Plummer küre modeli (astrofizik standardı)
    Disk,             // Yarı-düzenli: Disk şeklinde dağılım (galaksi modeli)
    Shell,            // Yarı-düzenli: Küresel kabuk
    DualGalaxy,       // Düzensiz: İki çarpışan galaksi
};

/**
 * Returns human-readable distribution name.
 */
const char* distribution_name(Distribution dist);

/**
 * Returns whether this distribution is uniform or irregular.
 */
const char* distribution_regularity(Distribution dist);

// ============================================================================
// Benchmark Configuration
// ============================================================================

struct BenchmarkConfig {
    // Particle counts to test
    std::vector<size_t> particle_counts = {1000, 5000, 10000, 50000, 100000};

    // Distributions to test
    std::vector<Distribution> distributions = {
        Distribution::Uniform,
        Distribution::Clustered,
        Distribution::Plummer,
        Distribution::Disk,
    };

    // Simulation parameters
    float theta = 0.5f;
    float dt = 0.001f;
    float epsilon_squared = 1e-10f;
    float gravity = 1.0f;
    int max_particles_per_leaf = 8;

    // Benchmark parameters (canonical protocol: --warmup 50 --steps 100)
    int warmup_steps = 10;      // Warmup iterations (not timed)
    int measured_steps = 30;    // Timed iterations, one sample each
    bool validate_correctness = true;  // Compare forces against a direct sum
    float max_l2_error = 2e-2f;        // Barnes-Hut: ||F - F_ref|| / ||F_ref||
    size_t validation_sample = 4096;   // Particles checked against the reference

    // Output
    std::string output_csv = "benchmark_results.csv";
    bool verbose = true;
};

// ============================================================================
// Benchmark Results
// ============================================================================

struct TimingResult {
    double mean_ms = 0;       // Mean execution time (milliseconds)
    double std_dev_ms = 0;    // Sample standard deviation
    double min_ms = 0;        // Minimum
    double max_ms = 0;        // Maximum
    double median_ms = 0;     // Median (primary statistic)
    double p25_ms = 0;        // First quartile
    double p75_ms = 0;        // Third quartile
    int num_samples = 0;
};

struct StepBreakdown {
    TimingResult tree_build;    // Tree construction (+ mass distribution)
    TimingResult mass_dist;     // Mass distribution (upward pass)
    TimingResult force_calc;    // Force calculation
    TimingResult integration;   // Particle integration
    TimingResult total;         // Total step time
};

struct CorrectnessResult {
    double max_relative_error = 0;    // Maximum per-particle relative force error
    double mean_relative_error = 0;   // Mean per-particle relative force error
    double l2_norm_error = 0;         // ||F - F_ref|| / ||F_ref|| (pass criterion)
    bool passed = true;               // Whether within acceptable tolerance
};

struct BenchmarkResult {
    // Identification
    std::string method_name;      // "CPU_Serial", "CPU_OpenMP", "GPU_Direct", "GPU_BH"
    Distribution distribution;
    size_t particle_count;

    // Timing
    StepBreakdown timing;

    // Correctness (compared to CPU serial reference)
    CorrectnessResult correctness;

    // Derived metrics
    double speedup_vs_serial = 0;     // Median-time speedup over CPU serial
    double gflops = 0;                // Not estimated (kept for CSV compatibility)
    double particles_per_second = 0;  // Throughput
};

// ============================================================================
// Data Generator
// ============================================================================

/**
 * Generates particle data for the specified distribution.
 * All arrays must be pre-allocated with at least `count` elements.
 */
void generate_distribution(
    Distribution dist,
    size_t count,
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    unsigned int seed = 42
);

// ============================================================================
// CPU Benchmark Functions
// ============================================================================

/**
 * Run CPU Barnes-Hut (octree) with one thread: warmup_steps untimed, then
 * num_steps timed kick-drift-kick steps. Positions, velocities and forces
 * of the final state are written back.
 */
StepBreakdown benchmark_cpu_serial(
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps,
    int warmup_steps = 0
);

/**
 * Same as benchmark_cpu_serial with all OpenMP threads.
 */
StepBreakdown benchmark_cpu_openmp(
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps,
    int warmup_steps = 0
);

/**
 * Time GPU O(N^2) force evaluations (the state is not advanced; the returned
 * forces belong to the input positions).
 */
StepBreakdown benchmark_gpu_direct(
    float* h_mass,
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps,
    int warmup_steps = 0
);

/**
 * Run GPU Barnes-Hut (LBVH) kick-drift-kick steps, optionally replayed from
 * a CUDA Graph. Final state is written back in the original order.
 */
StepBreakdown benchmark_gpu_bh(
    float* h_mass,
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps,
    int warmup_steps = 0,
    bool use_graph = false
);

/**
 * True when a CUDA device is usable.
 */
bool gpu_available();

/**
 * Barnes-Hut forces at the given state, without advancing time.
 * gpu_bh_forces returns false when no GPU is available.
 */
void cpu_bh_forces(const float* mass, const float* pos_x, const float* pos_y, const float* pos_z,
                   size_t count, const BenchmarkConfig& config,
                   float* force_x, float* force_y, float* force_z);
bool gpu_bh_forces(const float* mass, const float* pos_x, const float* pos_y, const float* pos_z,
                   size_t count, const BenchmarkConfig& config,
                   float* force_x, float* force_y, float* force_z);

// ============================================================================
// Correctness Validation
// ============================================================================

/**
 * Double-precision direct-sum forces for the particles listed in `indices`
 * (written to force_*[k] for indices[k]).
 */
void reference_forces_direct(
    const float* mass, const float* pos_x, const float* pos_y, const float* pos_z,
    size_t count, float epsilon_squared, float gravity,
    const std::vector<size_t>& indices,
    float* force_x, float* force_y, float* force_z
);

/**
 * Compare two force arrays and compute error metrics.
 */
CorrectnessResult validate_forces(
    const float* ref_fx, const float* ref_fy, const float* ref_fz,
    const float* test_fx, const float* test_fy, const float* test_fz,
    size_t count,
    float tolerance
);

// ============================================================================
// Benchmark Runner
// ============================================================================

/**
 * Run the complete benchmark suite.
 * Returns all results sorted by method and particle count.
 */
std::vector<BenchmarkResult> run_benchmark_suite(const BenchmarkConfig& config);

/**
 * Write results to CSV file for analysis.
 */
void write_results_csv(
    const std::vector<BenchmarkResult>& results,
    const std::string& filename
);

/**
 * Print formatted results table to stdout.
 */
void print_results_table(const std::vector<BenchmarkResult>& results);

/**
 * Print distribution-specific analysis.
 * Compares performance across regular vs irregular distributions.
 */
void print_distribution_analysis(const std::vector<BenchmarkResult>& results);

} // namespace benchmark
} // namespace barnes_hut
