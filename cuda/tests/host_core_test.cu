/**
 * Host-side verification of the GPU Barnes-Hut core (no GPU required).
 *
 * Executes the __host__ __device__ functions of bh_core.cuh in the same
 * order the kernels in bh_gpu.cu do, emulating thread concurrency where it
 * matters (refit leaves are processed in random order), and checks:
 *   1. split3_21 / morton63 against a bit-by-bit reference
 *   2. Karras hierarchy: valid binary tree, every leaf reached exactly once
 *   3. Refit: root mass/centre of mass, box containment for every node
 *   4. Forces against a double-precision direct sum (theta sweep,
 *      four distributions incl. coincident particles)
 *   5. Stack overflow is counted and the force stays finite
 *   6. Kick-drift-kick: bounded energy error on a two-body orbit
 *
 * Build: nvcc -std=c++20 -O2 -I cuda cuda/tests/host_core_test.cu
 */

#include "bh_core.cuh"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using namespace barnes_hut::gpu;

namespace {

int g_failures = 0;
int g_passes = 0;

void check(bool ok, const std::string& what) {
    std::printf("[%s] %s\n", ok ? "PASS" : "FAIL", what.c_str());
    (ok ? g_passes : g_failures)++;
}

// ---------------------------------------------------------------------------
// Particle sets
// ---------------------------------------------------------------------------

struct Particles {
    std::vector<float> x, y, z, m;
    explicit Particles(size_t n) : x(n), y(n), z(n), m(n) {}
    size_t size() const { return x.size(); }
};

Particles make_particles(const std::string& kind, size_t n, unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u(0.0f, 1.0f);
    std::normal_distribution<float> g(0.0f, 1.0f);
    Particles p(n);
    for (size_t i = 0; i < n; ++i) {
        p.m[i] = 0.5f + u(rng);
        if (kind == "uniform") {
            p.x[i] = 10.0f * u(rng); p.y[i] = 10.0f * u(rng); p.z[i] = 10.0f * u(rng);
        } else if (kind == "plummer") {
            // r = a / sqrt(X^{-2/3} - 1), a = 1, truncated at r = 20
            float r;
            do { r = 1.0f / std::sqrt(std::pow(std::max(u(rng), 1e-6f), -2.0f / 3.0f) - 1.0f); } while (r > 20.0f);
            const float ct = 2.0f * u(rng) - 1.0f, ph = 6.2831853f * u(rng);
            const float st = std::sqrt(1.0f - ct * ct);
            p.x[i] = r * st * std::cos(ph); p.y[i] = r * st * std::sin(ph); p.z[i] = r * ct;
        } else if (kind == "clustered") {
            const int c = static_cast<int>(i % 5);
            p.x[i] = 20.0f * c + 0.3f * g(rng); p.y[i] = 5.0f * (c % 2) + 0.3f * g(rng); p.z[i] = 0.3f * g(rng);
        } else {  // "coincident": 25% of the particles share four positions
            if (i % 4 == 0) {
                const float k = static_cast<float>(i % 16);
                p.x[i] = k; p.y[i] = 1.0f; p.z[i] = 2.0f;
            } else {
                p.x[i] = 16.0f * u(rng); p.y[i] = 16.0f * u(rng); p.z[i] = 16.0f * u(rng);
            }
        }
    }
    return p;
}

// ---------------------------------------------------------------------------
// Host emulation of the GPU pipeline
// ---------------------------------------------------------------------------

struct HostTree {
    int n = 0;
    std::vector<int> id;  // sorted index -> original index
    std::vector<float> x, y, z, m;
    std::vector<int> left, right, parent, visits;
    std::vector<float4> cm, lo, hi;

    TreeView view() {
        return TreeView{n, left.data(), right.data(), parent.data(), cm.data(), lo.data(), hi.data(), visits.data()};
    }
};

HostTree build(const Particles& p, float theta, unsigned shuffle_seed) {
    HostTree t;
    t.n = static_cast<int>(p.size());
    const int n = t.n;

    // bbox (cub::DeviceReduce over PointBox)
    Box box = empty_box();
    PointBox pb{p.x.data(), p.y.data(), p.z.data()};
    for (int i = 0; i < n; ++i) box = merge(box, pb(i));

    // keys + stable sort by key (cub::DeviceRadixSort is stable)
    std::vector<std::uint64_t> keys(n);
    for (int i = 0; i < n; ++i) keys[i] = morton_key(p.x[i], p.y[i], p.z[i], box);
    t.id.resize(n);
    std::iota(t.id.begin(), t.id.end(), 0);
    std::stable_sort(t.id.begin(), t.id.end(), [&](int a, int b) { return keys[a] < keys[b]; });
    std::vector<std::uint64_t> sorted(n);
    t.x.resize(n); t.y.resize(n); t.z.resize(n); t.m.resize(n);
    for (int k = 0; k < n; ++k) {
        const int s = t.id[k];
        sorted[k] = keys[s];
        t.x[k] = p.x[s]; t.y[k] = p.y[s]; t.z[k] = p.z[s]; t.m[k] = p.m[s];
    }

    // hierarchy
    t.left.assign(std::max(n - 1, 0), -1);
    t.right.assign(std::max(n - 1, 0), -1);
    t.parent.assign(2 * n - 1, -1);
    for (int i = 0; i < n - 1; ++i) build_internal_node(i, sorted.data(), n, t.left.data(), t.right.data(), t.parent.data());

    // refit in random leaf order
    t.visits.assign(std::max(n - 1, 0), 0);
    t.cm.assign(2 * n - 1, float4{0, 0, 0, 0});
    t.lo.assign(2 * n - 1, float4{0, 0, 0, 0});
    t.hi.assign(2 * n - 1, float4{0, 0, 0, 0});
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::shuffle(order.begin(), order.end(), std::mt19937(shuffle_seed));
    const float inv_theta2 = theta > 0.0f ? 1.0f / (theta * theta) : INFINITY;
    TreeView v = t.view();
    for (int j : order) refit_from_leaf(j, t.x.data(), t.y.data(), t.z.data(), t.m.data(), v, inv_theta2);
    return t;
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

void test_morton() {
    std::mt19937 rng(7);
    std::uniform_int_distribution<std::uint32_t> u(0, MORTON_MAX);
    bool ok = true;
    for (int k = 0; k < 200000 && ok; ++k) {
        const std::uint32_t a = u(rng), b = u(rng), c = u(rng);
        std::uint64_t ref = 0;
        for (int bit = 0; bit < 21; ++bit) {
            ref |= static_cast<std::uint64_t>((a >> bit) & 1u) << (3 * bit);
            ref |= static_cast<std::uint64_t>((b >> bit) & 1u) << (3 * bit + 1);
            ref |= static_cast<std::uint64_t>((c >> bit) & 1u) << (3 * bit + 2);
        }
        ok = (morton63(a, b, c) == ref);
    }
    check(ok, "morton63 == bit-by-bit reference (200k random triples)");
    check(quantize21(-1.0f, 0.0f, 1.0f) == 0 && quantize21(2.0f, 0.0f, 1.0f) == MORTON_MAX,
          "quantize21 clamps out-of-domain values");
}

int max_depth_and_check(const HostTree& t, bool& ok) {
    const int n = t.n;
    std::vector<int> seen(2 * n - 1, 0);
    int max_depth = 0;
    std::vector<std::pair<int, int>> stack{{0, 0}};
    while (!stack.empty()) {
        auto [node, depth] = stack.back();
        stack.pop_back();
        if (node < 0 || node >= 2 * n - 1) { ok = false; return -1; }
        seen[node]++;
        max_depth = std::max(max_depth, depth);
        if (node < n - 1) {
            if (t.parent[t.left[node]] != node || t.parent[t.right[node]] != node) ok = false;
            stack.push_back({t.left[node], depth + 1});
            stack.push_back({t.right[node], depth + 1});
        }
    }
    for (int k = 0; k < 2 * n - 1; ++k) ok = ok && seen[k] == 1;
    ok = ok && (n == 1 || t.parent[0] == -1);
    return max_depth;
}

void test_hierarchy_and_refit() {
    for (const char* kind : {"uniform", "plummer", "clustered", "coincident"}) {
        for (size_t n : {size_t(2), size_t(3), size_t(17), size_t(1000), size_t(50000)}) {
            const Particles p = make_particles(kind, n, 11);
            HostTree t = build(p, 0.5f, 3);
            bool ok = true;
            const int depth = max_depth_and_check(t, ok);

            // root mass / centre of mass vs double reference
            double M = 0, X = 0, Y = 0, Z = 0;
            for (size_t i = 0; i < n; ++i) { M += p.m[i]; X += p.m[i] * p.x[i]; Y += p.m[i] * p.y[i]; Z += p.m[i] * p.z[i]; }
            const float4 r = t.cm[0];
            const double scale = 1.0 + std::abs(X / M) + std::abs(Y / M) + std::abs(Z / M);
            const bool mass_ok = std::abs(r.w - M) <= 1e-5 * M &&
                                 std::abs(r.x - X / M) + std::abs(r.y - Y / M) + std::abs(r.z - Z / M) <= 1e-4 * scale;

            // every internal node box contains its children's boxes
            bool box_ok = true;
            for (int k = 0; k < t.n - 1; ++k) {
                for (int c : {t.left[k], t.right[k]}) {
                    box_ok = box_ok && t.lo[k].x <= t.lo[c].x && t.lo[k].y <= t.lo[c].y && t.lo[k].z <= t.lo[c].z &&
                             t.hi[k].x >= t.hi[c].x && t.hi[k].y >= t.hi[c].y && t.hi[k].z >= t.hi[c].z;
                }
            }
            check(ok && mass_ok && box_ok,
                  std::string("tree+refit ") + kind + " N=" + std::to_string(n) + " (depth " + std::to_string(depth) + ")");
        }
    }
}

struct ErrorStats {
    double median = 0, p99 = 0, max = 0;
    unsigned long long pp = 0, pc = 0, overflows = 0;
};

template <int STACK>
ErrorStats force_errors(const Particles& p, float theta, float eps2) {
    HostTree t = build(p, theta, 5);
    const int n = t.n;
    std::vector<double> fx(n), fy(n), fz(n);
    ErrorStats e;
    for (int k = 0; k < n; ++k) {
        TraversalCounts c;
        const float3 f = force_on<STACK>(k, t.left.data(), t.right.data(), t.cm.data(), t.lo.data(), t.hi.data(), n, eps2, 1.0f, c);
        const int o = t.id[k];
        fx[o] = f.x; fy[o] = f.y; fz[o] = f.z;
        e.pp += c.particle_particle; e.pc += c.particle_cell; e.overflows += c.stack_overflows;
    }
    std::vector<double> rel;
    for (int i = 0; i < n; ++i) {
        double rx = 0, ry = 0, rz = 0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            const double dx = double(p.x[j]) - p.x[i], dy = double(p.y[j]) - p.y[i], dz = double(p.z[j]) - p.z[i];
            const double r2 = dx * dx + dy * dy + dz * dz + eps2;
            const double s = double(p.m[i]) * p.m[j] / (r2 * std::sqrt(r2));
            rx += s * dx; ry += s * dy; rz += s * dz;
        }
        const double ref = std::sqrt(rx * rx + ry * ry + rz * rz);
        const double err = std::sqrt((fx[i] - rx) * (fx[i] - rx) + (fy[i] - ry) * (fy[i] - ry) + (fz[i] - rz) * (fz[i] - rz));
        rel.push_back(ref > 0 ? err / ref : err);
    }
    std::sort(rel.begin(), rel.end());
    e.median = rel[rel.size() / 2];
    e.p99 = rel[std::min(rel.size() - 1, rel.size() * 99 / 100)];
    e.max = rel.back();
    return e;
}

void test_forces() {
    const size_t n = 4000;
    std::printf("\n  %-11s %5s %10s %10s %10s %10s %10s\n", "dist", "theta", "median", "p99", "max", "P-P", "P-C");
    for (const char* kind : {"uniform", "plummer", "clustered", "coincident"}) {
        const Particles p = make_particles(kind, n, 21);
        // eps2 large enough that coincident pairs stay finite in float
        const float eps2 = 1e-4f;
        double previous_median = 0.0;
        for (float theta : {0.0f, 0.3f, 0.5f, 0.7f, 1.0f}) {
            const ErrorStats e = force_errors<64>(p, theta, eps2);
            std::printf("  %-11s %5.2f %10.3e %10.3e %10.3e %10llu %10llu\n", kind, theta, e.median, e.p99, e.max, e.pp, e.pc);
            if (theta == 0.0f) {
                // Coincident groups (b_max = 0) are exact monopoles even at theta = 0
                const bool exact_counts = std::string(kind) == "coincident" || (e.pp == n * (n - 1) && e.pc == 0);
                check(exact_counts && e.max < 1e-4,
                      std::string("theta=0 is an exact float direct sum (") + kind + ")");
            }
            if (theta == 0.5f) {
                // Regression guard; measured 4.8e-3..9.5e-3 / 4.1e-2..6.7e-2
                check(e.median < 1.2e-2 && e.p99 < 8e-2, std::string("theta=0.5 median<1.2e-2, p99<8e-2 (") + kind + ")");
            }
            check(e.median >= previous_median, std::string("median error grows with theta (") + kind + ")");
            previous_median = e.median;
            check(e.overflows == 0, std::string("no stack overflow, STACK=64 (") + kind + ", theta=" + std::to_string(theta) + ")");
        }
    }
}

void test_stack_overflow_guard() {
    const Particles p = make_particles("uniform", 2000, 3);
    const ErrorStats e = force_errors<4>(p, 0.0f, 1e-4f);
    check(e.overflows > 0 && std::isfinite(e.max), "STACK=4: overflows are counted, forces stay finite");
}

void test_leapfrog() {
    // Two equal masses on a circular orbit (G = 1, separation 1)
    const float v = std::sqrt(0.5f);
    float px[2] = {-0.5f, 0.5f}, py[2] = {0, 0}, pz[2] = {0, 0};
    float vx[2] = {0, 0}, vy[2] = {-v, v}, vz[2] = {0, 0};
    float fx[2], fy[2], fz[2], m[2] = {1, 1};
    const float eps2 = 0.0f, dt = 0.01f;

    auto forces = [&]() {
        Particles p(2);
        for (int i = 0; i < 2; ++i) { p.x[i] = px[i]; p.y[i] = py[i]; p.z[i] = pz[i]; p.m[i] = m[i]; }
        HostTree t = build(p, 0.5f, 1);
        for (int k = 0; k < 2; ++k) {
            TraversalCounts c;
            const float3 f = force_on<64>(k, t.left.data(), t.right.data(), t.cm.data(), t.lo.data(), t.hi.data(), 2, eps2, 1.0f, c);
            fx[t.id[k]] = f.x; fy[t.id[k]] = f.y; fz[t.id[k]] = f.z;
        }
    };
    auto energy = [&]() {
        const double dx = px[0] - px[1], dy = py[0] - py[1], dz = pz[0] - pz[1];
        double k = 0;
        for (int i = 0; i < 2; ++i) k += 0.5 * m[i] * (double(vx[i]) * vx[i] + double(vy[i]) * vy[i] + double(vz[i]) * vz[i]);
        return k - m[0] * m[1] / std::sqrt(dx * dx + dy * dy + dz * dz);
    };

    forces();
    const double e0 = energy();
    double worst = 0;
    for (int s = 0; s < 20000; ++s) {  // ~45 orbits
        for (int i = 0; i < 2; ++i) kick(i, vx, vy, vz, fx, fy, fz, m, 0.5f * dt);
        for (int i = 0; i < 2; ++i) drift(i, px, py, pz, vx, vy, vz, dt);
        forces();
        for (int i = 0; i < 2; ++i) kick(i, vx, vy, vz, fx, fy, fz, m, 0.5f * dt);
        worst = std::max(worst, std::abs((energy() - e0) / e0));
    }
    std::printf("  two-body KDK: max |dE/E| over 20000 steps = %.3e\n", worst);
    check(worst < 1e-4, "KDK leapfrog keeps |dE/E| < 1e-4 over ~45 orbits (float)");
}

} // namespace

int main() {
    std::printf("== Barnes-Hut GPU core: host emulation tests ==\n");
    test_morton();
    test_hierarchy_and_refit();
    test_forces();
    test_stack_overflow_guard();
    test_leapfrog();
    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures == 0 ? 0 : 1;
}
