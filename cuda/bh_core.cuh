#pragma once

/**
 * Barnes-Hut GPU core: pure __host__ __device__ building blocks.
 *
 * Every kernel in bh_gpu.cu is a thin grid-stride wrapper around one of the
 * functions below, so the whole algorithm (Morton keys, Karras hierarchy,
 * bottom-up refit, traversal, leapfrog) can be executed and verified on the
 * host without a GPU (see tests/host_core_test.cu).
 *
 * Tree layout (Karras 2012 linear BVH, binary radix tree over sorted keys):
 *   internal nodes  0 .. n-2      (node 0 is the root when n >= 2)
 *   leaf nodes      n-1 .. 2n-2   (leaf of sorted particle j is node n-1+j)
 * For every node:
 *   cm[node] = (center of mass xyz, total mass)
 *   lo[node] = (tight bounding box min xyz, (size / theta)^2)
 *   hi[node] = (tight bounding box max xyz, unused)
 */

#include <cuda/atomic>
#include <vector_types.h>
#include <cstdint>
#include <cfloat>
#include <cmath>

#if defined(__CUDACC__)
#define BH_HD __host__ __device__ __forceinline__
#else
#define BH_HD inline
#endif

namespace barnes_hut {
namespace gpu {

// ============================================================================
// Scalar helpers (identical results on host and device where it matters)
// ============================================================================

BH_HD int clz64(std::uint64_t x) {
#if defined(__CUDA_ARCH__)
    return __clzll(static_cast<long long>(x));
#else
    return x == 0 ? 64 : __builtin_clzll(x);
#endif
}

BH_HD int clz32(std::uint32_t x) {
#if defined(__CUDA_ARCH__)
    return __clz(static_cast<int>(x));
#else
    return x == 0 ? 32 : __builtin_clz(x);
#endif
}

BH_HD float rsqrt_approx(float x) {
#if defined(__CUDA_ARCH__)
    return rsqrtf(x);
#else
    return 1.0f / std::sqrt(x);
#endif
}

// ============================================================================
// Bounding box
// ============================================================================

struct Box {
    float3 lo;
    float3 hi;
};

BH_HD Box empty_box() {
    return Box{float3{FLT_MAX, FLT_MAX, FLT_MAX}, float3{-FLT_MAX, -FLT_MAX, -FLT_MAX}};
}

BH_HD Box merge(const Box& a, const Box& b) {
    return Box{float3{fminf(a.lo.x, b.lo.x), fminf(a.lo.y, b.lo.y), fminf(a.lo.z, b.lo.z)},
               float3{fmaxf(a.hi.x, b.hi.x), fmaxf(a.hi.y, b.hi.y), fmaxf(a.hi.z, b.hi.z)}};
}

struct BoxMerge {
    BH_HD Box operator()(const Box& a, const Box& b) const { return merge(a, b); }
};

// Maps a particle index to its degenerate box; fed to cub::DeviceReduce
// through cuda::transform_iterator, so no intermediate array is stored.
struct PointBox {
    const float* x;
    const float* y;
    const float* z;
    BH_HD Box operator()(int i) const {
        const float3 p{x[i], y[i], z[i]};
        return Box{p, p};
    }
};

// ============================================================================
// 63-bit Morton keys (21 bits per axis)
// ============================================================================

inline constexpr int MORTON_BITS_PER_AXIS = 21;
inline constexpr std::uint32_t MORTON_MAX = (1u << MORTON_BITS_PER_AXIS) - 1u;

// Inserts two zero bits after each of the low 21 bits: bit j -> bit 3j.
BH_HD std::uint64_t split3_21(std::uint32_t v) {
    std::uint64_t x = v & 0x1fffffu;
    x = (x | (x << 32)) & 0x001f00000000ffffull;
    x = (x | (x << 16)) & 0x001f0000ff0000ffull;
    x = (x | (x << 8))  & 0x100f00f00f00f00full;
    x = (x | (x << 4))  & 0x10c30c30c30c30c3ull;
    x = (x | (x << 2))  & 0x1249249249249249ull;
    return x;
}

BH_HD std::uint64_t morton63(std::uint32_t qx, std::uint32_t qy, std::uint32_t qz) {
    return split3_21(qx) | (split3_21(qy) << 1) | (split3_21(qz) << 2);
}

// Quantizes v into [0, 2^21 - 1] relative to a cubic domain; clamping
// happens in float before the cast, so the conversion is always defined.
BH_HD std::uint32_t quantize21(float v, float lo, float inv_extent) {
    const float t = fminf(fmaxf((v - lo) * inv_extent, 0.0f), 1.0f);
    return static_cast<std::uint32_t>(t * static_cast<float>(MORTON_MAX));
}

BH_HD std::uint64_t morton_key(float x, float y, float z, const Box& domain) {
    const float ex = domain.hi.x - domain.lo.x;
    const float ey = domain.hi.y - domain.lo.y;
    const float ez = domain.hi.z - domain.lo.z;
    const float extent = fmaxf(ex, fmaxf(ey, ez));
    const float inv = extent > 0.0f ? 1.0f / extent : 0.0f;
    return morton63(quantize21(x, domain.lo.x, inv),
                    quantize21(y, domain.lo.y, inv),
                    quantize21(z, domain.lo.z, inv));
}

// ============================================================================
// Karras hierarchy construction (one independent call per internal node)
// ============================================================================

// Length of the common prefix of keys i and j; equal keys are disambiguated
// by their (unique) sorted index, which keeps the tree well defined for
// coincident particles.
BH_HD int common_prefix(const std::uint64_t* keys, int n, int i, long long j) {
    if (j < 0 || j >= n) {
        return -1;
    }
    const std::uint64_t ki = keys[i];
    const std::uint64_t kj = keys[j];
    if (ki == kj) {
        return 64 + clz32(static_cast<std::uint32_t>(i) ^ static_cast<std::uint32_t>(j));
    }
    return clz64(ki ^ kj);
}

BH_HD void build_internal_node(int i, const std::uint64_t* keys, int n,
                               int* left, int* right, int* parent) {
    const int d = (common_prefix(keys, n, i, i + 1) - common_prefix(keys, n, i, i - 1)) >= 0 ? 1 : -1;

    // Upper bound for the range length, then binary search for the other end
    const int delta_min = common_prefix(keys, n, i, static_cast<long long>(i) - d);
    long long l_max = 2;
    while (common_prefix(keys, n, i, i + l_max * d) > delta_min) {
        l_max <<= 1;
    }
    long long l = 0;
    for (long long t = l_max >> 1; t >= 1; t >>= 1) {
        if (common_prefix(keys, n, i, i + (l + t) * d) > delta_min) {
            l += t;
        }
    }
    const long long j = i + l * d;

    // Split position: highest differing bit inside [i, j]
    const int delta_node = common_prefix(keys, n, i, j);
    long long s = 0;
    long long t = l;
    do {
        t = (t + 1) >> 1;
        if (common_prefix(keys, n, i, i + (s + t) * d) > delta_node) {
            s += t;
        }
    } while (t > 1);
    const int gamma = static_cast<int>(i + s * d + (d < 0 ? d : 0));

    const int range_lo = static_cast<int>(i < j ? i : j);
    const int range_hi = static_cast<int>(i < j ? j : i);
    const int a = (range_lo == gamma) ? (n - 1 + gamma) : gamma;
    const int b = (range_hi == gamma + 1) ? (n - 1 + gamma + 1) : (gamma + 1);

    left[i] = a;
    right[i] = b;
    parent[a] = i;
    parent[b] = i;
}

// ============================================================================
// Bottom-up refit (centre of mass, tight box, opening radius)
// ============================================================================

struct TreeView {
    int n;               // particle count; 2n-1 nodes
    const int* left;     // [n-1]
    const int* right;    // [n-1]
    const int* parent;   // [2n-1], -1 for the root
    float4* cm;          // [2n-1]
    float4* lo;          // [2n-1]
    float4* hi;          // [2n-1]
    int* visits;         // [n-1], zeroed before each refit
};

BH_HD void combine_children(int node, const TreeView& t, float inv_theta2) {
    const int a = t.left[node];
    const int b = t.right[node];
    const float4 ca = t.cm[a], cb = t.cm[b];
    const float4 la = t.lo[a], lb = t.lo[b];
    const float4 ha = t.hi[a], hb = t.hi[b];

    const float m = ca.w + cb.w;
    float4 c;
    if (m > 0.0f) {
        const float inv_m = 1.0f / m;
        c = float4{(ca.w * ca.x + cb.w * cb.x) * inv_m,
                   (ca.w * ca.y + cb.w * cb.y) * inv_m,
                   (ca.w * ca.z + cb.w * cb.z) * inv_m, m};
    } else {
        c = float4{0.5f * (ca.x + cb.x), 0.5f * (ca.y + cb.y), 0.5f * (ca.z + cb.z), 0.0f};
    }

    const float3 bl{fminf(la.x, lb.x), fminf(la.y, lb.y), fminf(la.z, lb.z)};
    const float3 bh{fmaxf(ha.x, hb.x), fmaxf(ha.y, hb.y), fmaxf(ha.z, hb.z)};

    const float size = fmaxf(bh.x - bl.x, fmaxf(bh.y - bl.y, bh.z - bl.z));

    // Accept the cell when |p - com| > size / theta (size = longest edge of
    // the tight box). At equal interaction count this was more accurate than
    // the b_max criterion (distance com -> farthest corner) for uniform,
    // Plummer and clustered sets in the host tests; the self-force case it
    // does not exclude is handled by the inside-box test in force_on().
    // theta = 0 gives an infinite radius (never accepted); size = 0
    // (coincident particles) is an exact monopole.
    const float open2 = size > 0.0f ? size * size * inv_theta2 : 0.0f;

    t.cm[node] = c;
    t.lo[node] = float4{bl.x, bl.y, bl.z, open2};
    t.hi[node] = float4{bh.x, bh.y, bh.z, 0.0f};
}

// Leaf j writes its own data and walks towards the root. The first thread
// that reaches an internal node stops; the second one combines both
// children. acq_rel ordering on the visit counter makes the sibling's node
// data (plain stores before its release) visible to the second thread.
BH_HD void refit_from_leaf(int j, const float* px, const float* py, const float* pz,
                           const float* mass, const TreeView& t, float inv_theta2) {
    const int leaf = t.n - 1 + j;
    t.cm[leaf] = float4{px[j], py[j], pz[j], mass[j]};
    t.lo[leaf] = float4{px[j], py[j], pz[j], 0.0f};
    t.hi[leaf] = float4{px[j], py[j], pz[j], 0.0f};

    int node = t.parent[leaf];
    while (node >= 0) {
        cuda::atomic_ref<int, cuda::thread_scope_device> visits(t.visits[node]);
        if (visits.fetch_add(1, cuda::std::memory_order_acq_rel) == 0) {
            return;
        }
        combine_children(node, t, inv_theta2);
        node = t.parent[node];
    }
}

// ============================================================================
// Force traversal
// ============================================================================

struct TraversalCounts {
    unsigned int particle_particle = 0;
    unsigned int particle_cell = 0;
    unsigned int stack_overflows = 0;
};

BH_HD void accumulate(float3& f, const float4& p, const float4& src, float eps2, float G) {
    const float dx = src.x - p.x;
    const float dy = src.y - p.y;
    const float dz = src.z - p.z;
    const float r2 = dx * dx + dy * dy + dz * dz + eps2;
    const float inv_r = rsqrt_approx(r2);
    const float s = G * p.w * src.w * (inv_r * inv_r * inv_r);
    f.x += s * dx;
    f.y += s * dy;
    f.z += s * dz;
}

BH_HD bool inside_box(const float4& p, const float4& lo, const float4& hi) {
    return p.x >= lo.x && p.x <= hi.x && p.y >= lo.y && p.y <= hi.y && p.z >= lo.z && p.z <= hi.z;
}

// Force on sorted particle i. Depth-first with an explicit stack: opening a
// node pushes its right child and continues with the left, so the stack
// holds at most one entry per tree level. If a level would exceed
// STACK_SIZE the node is used as a monopole instead (never silently dropped)
// and the event is counted, so the caller can detect it.
template <int STACK_SIZE>
BH_HD float3 force_on(int i, const int* left, const int* right, const float4* cm,
                      const float4* lo, const float4* hi, int n, float eps2, float G,
                      TraversalCounts& counts) {
    float3 f{0.0f, 0.0f, 0.0f};
    const int self = n - 1 + i;
    const float4 p = cm[self];
    if (n < 2) {
        return f;
    }

    int stack[STACK_SIZE];
    int top = 0;
    int node = 0;
    for (;;) {
        if (node >= n - 1) {
            if (node != self) {
                accumulate(f, p, cm[node], eps2, G);
                ++counts.particle_particle;
            }
        } else {
            const float4 c = cm[node];
            const float4 l = lo[node];
            const float dx = c.x - p.x, dy = c.y - p.y, dz = c.z - p.z;
            const float d2 = dx * dx + dy * dy + dz * dz;
            // A cell containing the particle is always opened (no self-force)
            if (d2 > l.w && !inside_box(p, l, hi[node])) {
                accumulate(f, p, c, eps2, G);
                ++counts.particle_cell;
            } else if (top < STACK_SIZE) {
                stack[top++] = right[node];
                node = left[node];
                continue;
            } else {
                accumulate(f, p, c, eps2, G);
                ++counts.particle_cell;
                ++counts.stack_overflows;
            }
        }
        if (top == 0) {
            break;
        }
        node = stack[--top];
    }
    return f;
}

// ============================================================================
// Kick-drift-kick leapfrog
// ============================================================================

BH_HD void kick(int i, float* vx, float* vy, float* vz, const float* fx, const float* fy,
                const float* fz, const float* mass, float h) {
    const float m = mass[i];
    if (m > 0.0f) {
        const float s = h / m;
        vx[i] += fx[i] * s;
        vy[i] += fy[i] * s;
        vz[i] += fz[i] * s;
    }
}

BH_HD void drift(int i, float* px, float* py, float* pz, const float* vx, const float* vy,
                 const float* vz, float dt) {
    px[i] += vx[i] * dt;
    py[i] += vy[i] * dt;
    pz[i] += vz[i] * dt;
}

} // namespace gpu
} // namespace barnes_hut
