#pragma once

/**
 * Benchmark Framework for Barnes-Hut Simulation
 *
 * Systematic comparison of:
 * 1. CPU Serial
 * 2. CPU OpenMP (parallel)
 * 3. GPU CUDA 13.2 (direct N-body)
 * 4. GPU CUDA 13.2 (Barnes-Hut tree)
 *
 * Features:
 * - Statistical rigor: multiple runs, warmup, std deviation
 * - Data distribution variants: uniform, clustered, Plummer, disk
 * - Scaling analysis: varying particle counts
 * - Correctness validation: force comparison (L2 norm)
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

    // Benchmark parameters
    int warmup_steps = 2;       // Warmup iterations (not timed)
    int measured_steps = 5;     // Timed iterations
    bool validate_correctness = true;  // Compare GPU results against CPU reference
    float max_relative_error = 1e-3f;  // Maximum acceptable relative error

    // Output
    std::string output_csv = "benchmark_results.csv";
    bool verbose = true;
};

// ============================================================================
// Benchmark Results
// ============================================================================

struct TimingResult {
    double mean_ms;       // Mean execution time (milliseconds)
    double std_dev_ms;    // Standard deviation
    double min_ms;        // Minimum
    double max_ms;        // Maximum
    double median_ms;     // Median
    int num_samples;
};

struct StepBreakdown {
    TimingResult bbox;          // Bounding box computation
    TimingResult tree_build;    // Tree construction
    TimingResult mass_dist;     // Mass distribution (upward pass)
    TimingResult force_calc;    // Force calculation
    TimingResult integration;   // Particle integration
    TimingResult total;         // Total step time
};

struct CorrectnessResult {
    double max_relative_error;    // Maximum relative force error
    double mean_relative_error;   // Mean relative force error
    double l2_norm_error;         // L2 norm of force difference
    bool passed;                  // Whether within acceptable tolerance
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
    double speedup_vs_serial;     // Speedup over CPU serial
    double gflops;                // Estimated GFLOP/s
    double particles_per_second;  // Throughput
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
 * Run CPU serial Barnes-Hut simulation.
 * Returns step breakdown timing.
 */
StepBreakdown benchmark_cpu_serial(
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps
);

/**
 * Run CPU OpenMP Barnes-Hut simulation.
 */
StepBreakdown benchmark_cpu_openmp(
    float* mass,
    float* pos_x, float* pos_y, float* pos_z,
    float* vel_x, float* vel_y, float* vel_z,
    float* force_x, float* force_y, float* force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps
);

/**
 * Run GPU direct N-body (O(N^2)) using CUDA 13.2.
 */
StepBreakdown benchmark_gpu_direct(
    float* h_mass,
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps
);

/**
 * Run GPU Barnes-Hut tree using CUDA 13.2 kernels.
 */
StepBreakdown benchmark_gpu_bh(
    float* h_mass,
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count,
    const BenchmarkConfig& config,
    int num_steps
);

// ============================================================================
// Correctness Validation
// ============================================================================

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
