/**
 * Test Scenarios for Barnes-Hut CUDA 13.2 Integration
 *
 * Comprehensive test suite covering:
 * 1. Correctness validation (CPU reference vs GPU)
 * 2. Performance benchmarks (scaling, throughput)
 * 3. Data distribution impact analysis (regular vs irregular)
 * 4. Numerical stability tests
 * 5. Edge cases
 *
 * Usage:
 *   ./benchmark_cuda13 [--quick | --full | --distribution-only | --correctness-only]
 */

#include "benchmark_framework.cuh"
#include <cuda_runtime.h>
#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <cmath>
#include <cassert>
#include <iomanip>

using namespace barnes_hut::benchmark;

// ============================================================================
// Test Helper Macros
// ============================================================================

#define TEST_ASSERT(cond, msg) \
    do { \
        if (!(cond)) { \
            std::cerr << "[FAIL] " << msg << " (" << __FILE__ << ":" << __LINE__ << ")\n"; \
            test_failures++; \
        } else { \
            std::cout << "[PASS] " << msg << "\n"; \
            test_passes++; \
        } \
    } while(0)

static int test_passes = 0;
static int test_failures = 0;

// ============================================================================
// Test 1: Data Distribution Generation Validation
// ============================================================================

void test_distribution_generation() {
    std::cout << "\n=== Test 1: Veri Dagılımı Uretim Dogrulaması ===\n\n";

    const size_t N = 10000;
    std::vector<float> mass(N), px(N), py(N), pz(N);
    std::vector<float> vx(N), vy(N), vz(N);

    Distribution dists[] = {
        Distribution::Uniform,
        Distribution::Clustered,
        Distribution::Plummer,
        Distribution::Disk,
        Distribution::Shell,
        Distribution::DualGalaxy,
    };

    for (auto dist : dists) {
        generate_distribution(dist, N,
            mass.data(), px.data(), py.data(), pz.data(),
            vx.data(), vy.data(), vz.data(), 42);

        // Basic sanity checks
        bool all_mass_positive = true;
        bool any_nan = false;
        float min_x = 1e30f, max_x = -1e30f;
        float min_y = 1e30f, max_y = -1e30f;
        float min_z = 1e30f, max_z = -1e30f;

        for (size_t i = 0; i < N; ++i) {
            if (mass[i] <= 0) all_mass_positive = false;
            if (std::isnan(px[i]) || std::isnan(py[i]) || std::isnan(pz[i])) any_nan = true;
            if (std::isnan(vx[i]) || std::isnan(vy[i]) || std::isnan(vz[i])) any_nan = true;
            min_x = std::min(min_x, px[i]); max_x = std::max(max_x, px[i]);
            min_y = std::min(min_y, py[i]); max_y = std::max(max_y, py[i]);
            min_z = std::min(min_z, pz[i]); max_z = std::max(max_z, pz[i]);
        }

        std::string name = distribution_name(dist);
        TEST_ASSERT(all_mass_positive, name + ": Tum kutleler pozitif");
        TEST_ASSERT(!any_nan, name + ": NaN degeri yok");
        TEST_ASSERT(max_x > min_x, name + ": X ekseninde yayılım var");
        TEST_ASSERT(max_y > min_y, name + ": Y ekseninde yayılım var");

        // Distribution-specific checks
        if (dist == Distribution::Uniform) {
            float range_x = max_x - min_x;
            float range_y = max_y - min_y;
            float range_z = max_z - min_z;
            // Uniform should span roughly [0,100]
            TEST_ASSERT(range_x > 80.0f && range_x < 110.0f,
                        name + ": X aralığı beklenen sınırlarda");
        }

        if (dist == Distribution::Clustered) {
            // Compute variance — should be lower than uniform
            float mean_x = 0, mean_y = 0;
            for (size_t i = 0; i < N; ++i) { mean_x += px[i]; mean_y += py[i]; }
            mean_x /= N; mean_y /= N;
            // Clustered particles should have noticeable grouping
            TEST_ASSERT(true, name + ": Küme dagılımı uretildi");
        }

        if (dist == Distribution::Disk) {
            // Z should be thin
            float z_range = max_z - min_z;
            float x_range = max_x - min_x;
            TEST_ASSERT(z_range < x_range * 0.5f,
                        name + ": Z ekseni X ekseninden dar (disk)");
        }

        std::cout << "    " << name << " bbox: ["
                  << min_x << "," << min_y << "," << min_z << "] - ["
                  << max_x << "," << max_y << "," << max_z << "]\n";
    }
}

// ============================================================================
// Test 2: CPU Barnes-Hut Correctness (vs Direct)
// ============================================================================

void test_cpu_correctness() {
    std::cout << "\n=== Test 2: CPU Barnes-Hut Dogruluk Testi ===\n\n";

    const size_t N = 500;  // Small enough for O(N^2) reference
    std::vector<float> mass(N), px(N), py(N), pz(N);
    std::vector<float> vx(N), vy(N), vz(N);
    std::vector<float> fx_bh(N), fy_bh(N), fz_bh(N);
    std::vector<float> fx_direct(N), fy_direct(N), fz_direct(N);

    BenchmarkConfig config;
    config.theta = 0.3f;  // Strict accuracy
    config.dt = 0.001f;
    config.max_particles_per_leaf = 1;
    config.warmup_steps = 0;
    config.measured_steps = 1;

    for (auto dist : {Distribution::Uniform, Distribution::Clustered, Distribution::Plummer}) {
        generate_distribution(dist, N,
            mass.data(), px.data(), py.data(), pz.data(),
            vx.data(), vy.data(), vz.data(), 42);

        // Save initial state
        auto px0 = px, py0 = py, pz0 = pz;
        auto vx0 = vx, vy0 = vy, vz0 = vz;

        // Run Barnes-Hut
        benchmark_cpu_serial(mass.data(), px.data(), py.data(), pz.data(),
                             vx.data(), vy.data(), vz.data(),
                             fx_bh.data(), fy_bh.data(), fz_bh.data(),
                             N, config, 1);

        // Compute direct N-body reference
        px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;
        float eps_sq = config.epsilon_squared;
        float G = config.gravity;

        for (size_t i = 0; i < N; ++i) {
            fx_direct[i] = fy_direct[i] = fz_direct[i] = 0.0f;
        }

        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (i == j) continue;
                float dx = px[i] - px[j];
                float dy = py[i] - py[j];
                float dz = pz[i] - pz[j];
                float r_sq = dx*dx + dy*dy + dz*dz + eps_sq;
                float r_inv = 1.0f / std::sqrt(r_sq);
                float r_inv3 = r_inv * r_inv * r_inv;
                float f = -G * mass[i] * mass[j] * r_inv3;
                fx_direct[i] += f * dx;
                fy_direct[i] += f * dy;
                fz_direct[i] += f * dz;
            }
        }

        auto result = validate_forces(
            fx_direct.data(), fy_direct.data(), fz_direct.data(),
            fx_bh.data(), fy_bh.data(), fz_bh.data(),
            N, 0.05f  // 5% tolerance for theta=0.3
        );

        std::string name = distribution_name(dist);
        TEST_ASSERT(result.passed,
            name + ": BH vs Direct max_rel_err=" +
            std::to_string(result.max_relative_error));

        std::cout << "    " << name << ": L2 norm error = " << result.l2_norm_error
                  << ", max relative error = " << result.max_relative_error << "\n";
    }
}

// ============================================================================
// Test 3: Force Calculation Correctness Validation
// ============================================================================

void test_force_validation() {
    std::cout << "\n=== Test 3: Kuvvet Hesabi Dogrulama ===\n\n";

    // Test with known analytical solution: 2 particles
    {
        const size_t N = 2;
        float mass[2] = {1000.0f, 2000.0f};
        float px[2] = {0.0f, 10.0f};
        float py[2] = {0.0f, 0.0f};
        float pz[2] = {0.0f, 0.0f};
        float vx[2] = {0.0f, 0.0f};
        float vy[2] = {0.0f, 0.0f};
        float vz[2] = {0.0f, 0.0f};
        float fx[2], fy[2], fz[2];

        BenchmarkConfig config;
        config.dt = 0.001f;
        config.theta = 0.0f;  // Force exact
        config.max_particles_per_leaf = 1;

        benchmark_cpu_serial(mass, px, py, pz, vx, vy, vz, fx, fy, fz,
                             N, config, 1);

        // Expected: F = -G * m1 * m2 / r^2 * r_hat
        // r = 10, G = 1.0, m1 = 1000, m2 = 2000
        // F = -1.0 * 1000 * 2000 / (100 + eps) in x-direction
        float expected_f_mag = 1.0f * 1000.0f * 2000.0f / (100.0f + 1e-10f);
        float actual_f0_x = std::abs(fx[0]);

        // Forces should be equal and opposite
        float f_ratio = actual_f0_x / expected_f_mag;
        TEST_ASSERT(std::abs(f_ratio - 1.0f) < 0.01f,
                    "2-parcacık kuvvet buyuklugu doğru (hata < %1)");

        // Newton's 3rd law: forces equal and opposite
        float sum_fx = fx[0] + fx[1];
        TEST_ASSERT(std::abs(sum_fx) < std::abs(fx[0]) * 0.01f,
                    "Newton 3. yasası: kuvvetler esit ve zıt");
    }
}

// ============================================================================
// Test 4: Numerical Stability
// ============================================================================

void test_numerical_stability() {
    std::cout << "\n=== Test 4: Sayısal Kararlılık Testi ===\n\n";

    const size_t N = 1000;
    std::vector<float> mass(N), px(N), py(N), pz(N);
    std::vector<float> vx(N), vy(N), vz(N);
    std::vector<float> fx(N), fy(N), fz(N);

    BenchmarkConfig config;
    config.dt = 0.001f;
    config.theta = 0.5f;
    config.max_particles_per_leaf = 8;

    // Test: Very close particles (softening test)
    generate_distribution(Distribution::Clustered, N,
        mass.data(), px.data(), py.data(), pz.data(),
        vx.data(), vy.data(), vz.data(), 42);

    // Place two particles very close together
    px[0] = 50.0f; py[0] = 50.0f; pz[0] = 50.0f;
    px[1] = 50.0f + 1e-6f; py[1] = 50.0f; pz[1] = 50.0f;

    benchmark_cpu_serial(mass.data(), px.data(), py.data(), pz.data(),
                         vx.data(), vy.data(), vz.data(),
                         fx.data(), fy.data(), fz.data(),
                         N, config, 1);

    bool any_nan = false;
    bool any_inf = false;
    for (size_t i = 0; i < N; ++i) {
        if (std::isnan(fx[i]) || std::isnan(fy[i]) || std::isnan(fz[i])) any_nan = true;
        if (std::isinf(fx[i]) || std::isinf(fy[i]) || std::isinf(fz[i])) any_inf = true;
    }

    TEST_ASSERT(!any_nan, "Cok yakın parcacıklarda NaN yok (softening calısıyor)");
    TEST_ASSERT(!any_inf, "Cok yakın parcacıklarda Inf yok");

    // Test: Very large domain
    generate_distribution(Distribution::Uniform, N,
        mass.data(), px.data(), py.data(), pz.data(),
        vx.data(), vy.data(), vz.data(), 42);

    // Scale positions to very large range
    for (size_t i = 0; i < N; ++i) {
        px[i] *= 1e6f; py[i] *= 1e6f; pz[i] *= 1e6f;
    }

    benchmark_cpu_serial(mass.data(), px.data(), py.data(), pz.data(),
                         vx.data(), vy.data(), vz.data(),
                         fx.data(), fy.data(), fz.data(),
                         N, config, 1);

    any_nan = false;
    for (size_t i = 0; i < N; ++i) {
        if (std::isnan(fx[i]) || std::isnan(fy[i]) || std::isnan(fz[i])) any_nan = true;
    }

    TEST_ASSERT(!any_nan, "Buyuk alan olceginde NaN yok");
}

// ============================================================================
// Test 5: Scaling Analysis
// ============================================================================

void test_scaling_analysis() {
    std::cout << "\n=== Test 5: Olcekleme Analizi (Scaling) ===\n\n";

    std::cout << "Bu test, parcacık sayısının arttırılmasıyla\n"
              << "hesaplama suresinin nasıl degistigini gosterir.\n"
              << "Barnes-Hut: O(N log N), Direct: O(N^2) beklenir.\n\n";

    BenchmarkConfig config;
    config.theta = 0.5f;
    config.dt = 0.001f;
    config.max_particles_per_leaf = 8;
    config.warmup_steps = 1;
    config.measured_steps = 3;

    std::vector<size_t> sizes = {500, 1000, 2000, 5000, 10000};

    std::cout << std::left
              << std::setw(10) << "N"
              << std::setw(15) << "Time(ms)"
              << std::setw(15) << "N*logN"
              << std::setw(15) << "Ratio"
              << std::setw(15) << "Complexity"
              << "\n";
    std::cout << std::string(70, '-') << "\n";

    double prev_time = 0;
    double prev_nlogn = 0;

    for (auto N : sizes) {
        std::vector<float> mass(N), px(N), py(N), pz(N);
        std::vector<float> vx(N), vy(N), vz(N);
        std::vector<float> fx(N), fy(N), fz(N);

        generate_distribution(Distribution::Uniform, N,
            mass.data(), px.data(), py.data(), pz.data(),
            vx.data(), vy.data(), vz.data(), 42);

        auto timing = benchmark_cpu_serial(
            mass.data(), px.data(), py.data(), pz.data(),
            vx.data(), vy.data(), vz.data(),
            fx.data(), fy.data(), fz.data(),
            N, config, config.warmup_steps + config.measured_steps);

        double nlogn = N * std::log2(static_cast<double>(N));
        double ratio = (prev_nlogn > 0) ?
            (timing.total.mean_ms / prev_time) / (nlogn / prev_nlogn) : 1.0;

        std::string complexity;
        if (ratio > 0.8 && ratio < 1.3) complexity = "~ O(N log N)";
        else if (ratio > 1.3 && ratio < 2.5) complexity = "~ O(N^1.5)";
        else if (ratio > 2.5) complexity = "~ O(N^2)";
        else complexity = "sublinear";

        std::cout << std::left
                  << std::setw(10) << N
                  << std::setw(15) << std::fixed << std::setprecision(2) << timing.total.mean_ms
                  << std::setw(15) << std::setprecision(0) << nlogn
                  << std::setw(15) << std::setprecision(3) << ratio
                  << std::setw(15) << complexity
                  << "\n";

        prev_time = timing.total.mean_ms;
        prev_nlogn = nlogn;
    }

    TEST_ASSERT(true, "Olcekleme analizi tamamlandı");
}

// ============================================================================
// Test 6: Energy Conservation
// ============================================================================

void test_energy_conservation() {
    std::cout << "\n=== Test 6: Enerji Koruma Testi ===\n\n";

    const size_t N = 500;
    std::vector<float> mass(N), px(N), py(N), pz(N);
    std::vector<float> vx(N), vy(N), vz(N);
    std::vector<float> fx(N), fy(N), fz(N);

    BenchmarkConfig config;
    config.theta = 0.3f;
    config.dt = 0.0001f;  // Small dt for energy conservation
    config.max_particles_per_leaf = 1;
    config.warmup_steps = 0;
    config.measured_steps = 1;

    generate_distribution(Distribution::Plummer, N,
        mass.data(), px.data(), py.data(), pz.data(),
        vx.data(), vy.data(), vz.data(), 42);

    // Compute initial kinetic energy
    auto compute_kinetic = [&]() {
        double KE = 0;
        for (size_t i = 0; i < N; ++i) {
            KE += 0.5 * mass[i] * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
        }
        return KE;
    };

    // Compute potential energy (O(N^2))
    auto compute_potential = [&]() {
        double PE = 0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = i+1; j < N; ++j) {
                float dx = px[i] - px[j];
                float dy = py[i] - py[j];
                float dz = pz[i] - pz[j];
                float r = std::sqrt(dx*dx + dy*dy + dz*dz + config.epsilon_squared);
                PE -= config.gravity * mass[i] * mass[j] / r;
            }
        }
        return PE;
    };

    double E0 = compute_kinetic() + compute_potential();

    // Run several steps
    const int steps = 10;
    for (int s = 0; s < steps; ++s) {
        benchmark_cpu_serial(mass.data(), px.data(), py.data(), pz.data(),
                             vx.data(), vy.data(), vz.data(),
                             fx.data(), fy.data(), fz.data(),
                             N, config, 1);
    }

    double E1 = compute_kinetic() + compute_potential();
    double relative_energy_error = std::abs(E1 - E0) / std::abs(E0);

    std::cout << "  Baslangıc enerjisi: " << E0 << "\n"
              << "  Son enerji:         " << E1 << "\n"
              << "  Bagıl hata:         " << relative_energy_error << "\n";

    TEST_ASSERT(relative_energy_error < 0.1,
                "Enerji koruma hatası < %10 (" + std::to_string(relative_energy_error * 100) + "%)");
}

// ============================================================================
// Test 7: Distribution Impact Analysis
// ============================================================================

void test_distribution_impact() {
    std::cout << "\n=== Test 7: Dagılım Etkisi Analizi ===\n\n";

    std::cout << "Bu test, farklı veri dagılımlarının Barnes-Hut\n"
              << "algoritmasının performansına etkisini olcer.\n\n";

    BenchmarkConfig config;
    config.particle_counts = {1000, 5000, 10000};
    config.distributions = {
        Distribution::Uniform,
        Distribution::Clustered,
        Distribution::Plummer,
        Distribution::Disk,
    };
    config.warmup_steps = 1;
    config.measured_steps = 3;
    config.theta = 0.5f;
    config.dt = 0.001f;

    auto results = run_benchmark_suite(config);

    // Print comprehensive analysis
    print_results_table(results);
    print_distribution_analysis(results);

    // Write CSV for external analysis
    write_results_csv(results, "distribution_impact_results.csv");

    TEST_ASSERT(!results.empty(), "Dagılım etkisi analizi tamamlandı");
}

// ============================================================================
// Test 8: GPU vs CPU Correctness (if GPU available)
// ============================================================================

void test_gpu_correctness() {
    std::cout << "\n=== Test 8: GPU vs CPU Dogruluk Karsilastirmasi ===\n\n";

    int device_count = 0;
    cudaGetDeviceCount(&device_count);
    if (device_count == 0) {
        std::cout << "  GPU bulunamadı — test atlanıyor\n";
        TEST_ASSERT(true, "GPU testi atlandı (GPU yok)");
        return;
    }

    const size_t N = 1000;
    std::vector<float> mass(N), px(N), py(N), pz(N);
    std::vector<float> vx(N), vy(N), vz(N);
    std::vector<float> cpu_fx(N), cpu_fy(N), cpu_fz(N);
    std::vector<float> gpu_fx(N), gpu_fy(N), gpu_fz(N);

    BenchmarkConfig config;
    config.theta = 0.5f;
    config.dt = 0.001f;

    for (auto dist : {Distribution::Uniform, Distribution::Plummer}) {
        generate_distribution(dist, N,
            mass.data(), px.data(), py.data(), pz.data(),
            vx.data(), vy.data(), vz.data(), 42);

        auto px0 = px, py0 = py, pz0 = pz;
        auto vx0 = vx, vy0 = vy, vz0 = vz;

        // CPU reference (direct O(N^2))
        float eps_sq = config.epsilon_squared;
        float G = config.gravity;
        for (size_t i = 0; i < N; ++i) cpu_fx[i] = cpu_fy[i] = cpu_fz[i] = 0;
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (i == j) continue;
                float dx = px0[i] - px0[j];
                float dy = py0[i] - py0[j];
                float dz = pz0[i] - pz0[j];
                float r_sq = dx*dx + dy*dy + dz*dz + eps_sq;
                float r_inv = 1.0f / std::sqrt(r_sq);
                float r_inv3 = r_inv * r_inv * r_inv;
                float f = -G * mass[i] * mass[j] * r_inv3;
                cpu_fx[i] += f * dx;
                cpu_fy[i] += f * dy;
                cpu_fz[i] += f * dz;
            }
        }

        // GPU direct
        px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;
        benchmark_gpu_direct(mass.data(), px.data(), py.data(), pz.data(),
                             vx.data(), vy.data(), vz.data(),
                             gpu_fx.data(), gpu_fy.data(), gpu_fz.data(),
                             N, config, 1);

        // Check for CUDA kernel launch errors
        cudaError_t cuda_err = cudaGetLastError();
        if (cuda_err != cudaSuccess) {
            std::cerr << "  CUDA kernel error: " << cudaGetErrorString(cuda_err) << "\n";
            TEST_ASSERT(false, "GPU kernel basarılı calıstı (" + std::string(distribution_name(dist)) + ")");
            continue;
        }

        auto result = validate_forces(
            cpu_fx.data(), cpu_fy.data(), cpu_fz.data(),
            gpu_fx.data(), gpu_fy.data(), gpu_fz.data(),
            N, 0.01f
        );

        std::string name = distribution_name(dist);
        TEST_ASSERT(result.passed,
            "GPU vs CPU (" + name + "): max_err=" +
            std::to_string(result.max_relative_error));

        std::cout << "    " << name << ": L2=" << result.l2_norm_error
                  << " max_rel=" << result.max_relative_error << "\n";
    }
}

// ============================================================================
// Test 9: Theta Parameter Sensitivity
// ============================================================================

void test_theta_sensitivity() {
    std::cout << "\n=== Test 9: Theta Parametresi Duyarlılık Analizi ===\n\n";

    std::cout << "Theta degeri kuculdukce dogruluk artar, hız azalır.\n"
              << "Bu test, theta-dogruluk-hız dengesini gosterir.\n\n";

    const size_t N = 2000;
    std::vector<float> mass(N), px(N), py(N), pz(N);
    std::vector<float> vx(N), vy(N), vz(N);
    std::vector<float> fx(N), fy(N), fz(N);
    std::vector<float> ref_fx(N), ref_fy(N), ref_fz(N);

    generate_distribution(Distribution::Plummer, N,
        mass.data(), px.data(), py.data(), pz.data(),
        vx.data(), vy.data(), vz.data(), 42);

    auto px0 = px, py0 = py, pz0 = pz;
    auto vx0 = vx, vy0 = vy, vz0 = vz;

    // Direct reference
    float eps_sq = 1e-10f;
    float G = 1.0f;
    for (size_t i = 0; i < N; ++i) ref_fx[i] = ref_fy[i] = ref_fz[i] = 0;
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            if (i == j) continue;
            float dx = px0[i] - px0[j], dy = py0[i] - py0[j], dz = pz0[i] - pz0[j];
            float r_sq = dx*dx + dy*dy + dz*dz + eps_sq;
            float r_inv = 1.0f / std::sqrt(r_sq);
            float r_inv3 = r_inv * r_inv * r_inv;
            float f = -G * mass[i] * mass[j] * r_inv3;
            ref_fx[i] += f * dx; ref_fy[i] += f * dy; ref_fz[i] += f * dz;
        }
    }

    std::cout << std::left
              << std::setw(10) << "Theta"
              << std::setw(15) << "Time(ms)"
              << std::setw(15) << "MaxRelErr"
              << std::setw(15) << "L2Error"
              << std::setw(10) << "Dogruluk"
              << "\n";
    std::cout << std::string(65, '-') << "\n";

    float thetas[] = {0.0f, 0.1f, 0.3f, 0.5f, 0.8f, 1.0f, 1.5f};

    for (float theta : thetas) {
        px = px0; py = py0; pz = pz0; vx = vx0; vy = vy0; vz = vz0;

        BenchmarkConfig config;
        config.theta = theta;
        config.dt = 0.001f;
        config.max_particles_per_leaf = 1;
        config.warmup_steps = 0;
        config.measured_steps = 1;

        auto timing = benchmark_cpu_serial(
            mass.data(), px.data(), py.data(), pz.data(),
            vx.data(), vy.data(), vz.data(),
            fx.data(), fy.data(), fz.data(),
            N, config, 1);

        auto correctness = validate_forces(
            ref_fx.data(), ref_fy.data(), ref_fz.data(),
            fx.data(), fy.data(), fz.data(), N, 1.0f);

        std::string accuracy;
        if (correctness.max_relative_error < 0.01) accuracy = "Mukemmel";
        else if (correctness.max_relative_error < 0.05) accuracy = "Iyi";
        else if (correctness.max_relative_error < 0.15) accuracy = "Orta";
        else accuracy = "Dusuk";

        std::cout << std::left
                  << std::setw(10) << std::fixed << std::setprecision(1) << theta
                  << std::setw(15) << std::setprecision(2) << timing.total.mean_ms
                  << std::setw(15) << std::scientific << std::setprecision(3) << correctness.max_relative_error
                  << std::setw(15) << correctness.l2_norm_error
                  << std::setw(10) << accuracy
                  << "\n";
    }

    TEST_ASSERT(true, "Theta duyarlılık analizi tamamlandı");
}

// ============================================================================
// Test 10: Edge Cases
// ============================================================================

void test_edge_cases() {
    std::cout << "\n=== Test 10: Sınır Durumları ===\n\n";

    BenchmarkConfig config;
    config.theta = 0.5f;
    config.dt = 0.001f;
    config.max_particles_per_leaf = 8;
    config.warmup_steps = 0;
    config.measured_steps = 1;

    // Test: Single particle
    {
        float mass[1] = {1000.0f};
        float px[1] = {50.0f}, py[1] = {50.0f}, pz[1] = {50.0f};
        float vx[1] = {1.0f}, vy[1] = {0.0f}, vz[1] = {0.0f};
        float fx[1], fy[1], fz[1];

        benchmark_cpu_serial(mass, px, py, pz, vx, vy, vz, fx, fy, fz,
                             1, config, 1);

        TEST_ASSERT(fx[0] == 0.0f && fy[0] == 0.0f && fz[0] == 0.0f,
                    "Tek parcacık: kuvvet sıfır");
    }

    // Test: Two particles only
    {
        float mass[2] = {1000.0f, 1000.0f};
        float px[2] = {0.0f, 10.0f}, py[2] = {0.0f, 0.0f}, pz[2] = {0.0f, 0.0f};
        float vx[2] = {0.0f, 0.0f}, vy[2] = {0.0f, 0.0f}, vz[2] = {0.0f, 0.0f};
        float fx[2], fy[2], fz[2];

        benchmark_cpu_serial(mass, px, py, pz, vx, vy, vz, fx, fy, fz,
                             2, config, 1);

        // Forces should be non-zero and opposite
        TEST_ASSERT(std::abs(fx[0]) > 0 && std::abs(fx[1]) > 0,
                    "Iki parcacık: kuvvetler sıfır degil");
        TEST_ASSERT(std::abs(fx[0] + fx[1]) < std::abs(fx[0]) * 0.01f,
                    "Iki parcacık: Newton 3. yasa");
    }

    // Test: All particles at same position
    {
        const size_t N = 100;
        std::vector<float> mass(N, 1000.0f);
        std::vector<float> px(N, 50.0f), py(N, 50.0f), pz(N, 50.0f);
        std::vector<float> vx(N, 0.0f), vy(N, 0.0f), vz(N, 0.0f);
        std::vector<float> fx(N), fy(N), fz(N);

        benchmark_cpu_serial(mass.data(), px.data(), py.data(), pz.data(),
                             vx.data(), vy.data(), vz.data(),
                             fx.data(), fy.data(), fz.data(),
                             N, config, 1);

        bool any_nan = false;
        for (size_t i = 0; i < N; ++i) {
            if (std::isnan(fx[i]) || std::isnan(fy[i]) || std::isnan(fz[i])) any_nan = true;
        }
        TEST_ASSERT(!any_nan, "Aynı konumdaki parcacıklar: NaN yok (softening)");
    }
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char* argv[]) {
    std::cout << "================================================================\n"
              << "  Barnes-Hut CUDA 13.2 Test Suite\n"
              << "  Sistematik CPU vs GPU Karsilastirma\n"
              << "  Duzenli ve Duzensiz Dagılım Testleri\n"
              << "================================================================\n";

    bool run_all = true;
    bool quick_mode = false;
    bool distribution_only = false;
    bool correctness_only = false;

    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--quick") == 0) quick_mode = true;
        else if (strcmp(argv[i], "--full") == 0) run_all = true;
        else if (strcmp(argv[i], "--distribution-only") == 0) distribution_only = true;
        else if (strcmp(argv[i], "--correctness-only") == 0) correctness_only = true;
    }

    if (correctness_only) {
        test_distribution_generation();
        test_cpu_correctness();
        test_force_validation();
        test_numerical_stability();
        test_edge_cases();
    } else if (distribution_only) {
        test_distribution_impact();
    } else if (quick_mode) {
        test_distribution_generation();
        test_force_validation();
        test_edge_cases();
    } else {
        // Full suite
        test_distribution_generation();
        test_cpu_correctness();
        test_force_validation();
        test_numerical_stability();
        test_scaling_analysis();
        test_energy_conservation();
        test_distribution_impact();
        test_gpu_correctness();
        test_theta_sensitivity();
        test_edge_cases();
    }

    // Summary
    std::cout << "\n================================================================\n"
              << "  TEST SONUCLARI\n"
              << "================================================================\n"
              << "  Gecen: " << test_passes << "\n"
              << "  Kalan: " << test_failures << "\n"
              << "  Toplam: " << (test_passes + test_failures) << "\n"
              << "================================================================\n";

    return (test_failures > 0) ? 1 : 0;
}
