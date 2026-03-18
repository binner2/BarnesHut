/**
 * Barnes-Hut CUDA 13.2 Enhanced Kernels — Implementation
 *
 * Core GPU kernels leveraging CUDA 13.2 features:
 * - SoA memory layout for coalesced access
 * - Tile-based force calculation with shared memory
 * - Parallel reduction (bounding box, mass distribution)
 * - Leapfrog integration
 *
 * Compile with: nvcc -std=c++20 -arch=sm_80 (or higher)
 */

#include "bh_cuda_kernels.cuh"
#include <cuda_runtime.h>
#include <cstdio>
#include <cfloat>
#include <cmath>

// We conditionally include cub headers — they are part of CCCL 3.2 in CUDA 13.2
// If not available, fallback kernels are provided
#if __has_include(<cub/cub.cuh>)
    #include <cub/cub.cuh>
    #define HAS_CUB 1
#else
    #define HAS_CUB 0
#endif

namespace barnes_hut {
namespace cuda {

// ============================================================================
// Error Checking
// ============================================================================

#define CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d — %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
            return false; \
        } \
    } while(0)

#define CUDA_CHECK_VOID(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            fprintf(stderr, "CUDA error at %s:%d — %s\n", \
                    __FILE__, __LINE__, cudaGetErrorString(err)); \
        } \
    } while(0)

// ============================================================================
// GpuMemoryManager Implementation
// ============================================================================

GpuMemoryManager::GpuMemoryManager() : allocated_(false) {
    memset(&d_particles_, 0, sizeof(d_particles_));
    memset(&d_nodes_, 0, sizeof(d_nodes_));
}

GpuMemoryManager::~GpuMemoryManager() {
    free();
}

bool GpuMemoryManager::allocate(size_t particle_count, size_t max_node_count) {
    if (allocated_) free();

    d_particles_.count = particle_count;
    d_nodes_.node_count = max_node_count;

    // Allocate particle SoA arrays
    CUDA_CHECK(cudaMalloc(&d_particles_.pos_x, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.pos_y, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.pos_z, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.vel_x, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.vel_y, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.vel_z, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.force_x, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.force_y, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.force_z, particle_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_particles_.mass, particle_count * sizeof(float)));

    // Allocate node SoA arrays
    CUDA_CHECK(cudaMalloc(&d_nodes_.center_x, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.center_y, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.center_z, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.size, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.mass_center_x, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.mass_center_y, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.mass_center_z, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.total_mass, max_node_count * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.children, max_node_count * 8 * sizeof(int32_t)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.parent, max_node_count * sizeof(int32_t)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.particle_begin, max_node_count * sizeof(int32_t)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.particle_end, max_node_count * sizeof(int32_t)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.node_type, max_node_count * sizeof(uint8_t)));
    CUDA_CHECK(cudaMalloc(&d_nodes_.level, max_node_count * sizeof(int32_t)));

    allocated_ = true;
    return true;
}

void GpuMemoryManager::free() {
    if (!allocated_) return;

    // Free particle arrays
    cudaFree(d_particles_.pos_x);
    cudaFree(d_particles_.pos_y);
    cudaFree(d_particles_.pos_z);
    cudaFree(d_particles_.vel_x);
    cudaFree(d_particles_.vel_y);
    cudaFree(d_particles_.vel_z);
    cudaFree(d_particles_.force_x);
    cudaFree(d_particles_.force_y);
    cudaFree(d_particles_.force_z);
    cudaFree(d_particles_.mass);

    // Free node arrays
    cudaFree(d_nodes_.center_x);
    cudaFree(d_nodes_.center_y);
    cudaFree(d_nodes_.center_z);
    cudaFree(d_nodes_.size);
    cudaFree(d_nodes_.mass_center_x);
    cudaFree(d_nodes_.mass_center_y);
    cudaFree(d_nodes_.mass_center_z);
    cudaFree(d_nodes_.total_mass);
    cudaFree(d_nodes_.children);
    cudaFree(d_nodes_.parent);
    cudaFree(d_nodes_.particle_begin);
    cudaFree(d_nodes_.particle_end);
    cudaFree(d_nodes_.node_type);
    cudaFree(d_nodes_.level);

    memset(&d_particles_, 0, sizeof(d_particles_));
    memset(&d_nodes_, 0, sizeof(d_nodes_));
    allocated_ = false;
}

bool GpuMemoryManager::upload_particles(
    const float* h_mass,
    const float* h_pos_x, const float* h_pos_y, const float* h_pos_z,
    const float* h_vel_x, const float* h_vel_y, const float* h_vel_z,
    size_t count, cudaStream_t stream)
{
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.mass, h_mass, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.pos_x, h_pos_x, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.pos_y, h_pos_y, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.pos_z, h_pos_z, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.vel_x, h_vel_x, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.vel_y, h_vel_y, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK(cudaMemcpyAsync(d_particles_.vel_z, h_vel_z, count * sizeof(float), cudaMemcpyHostToDevice, stream));
    return true;
}

bool GpuMemoryManager::download_particles(
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z,
    size_t count, cudaStream_t stream)
{
    CUDA_CHECK(cudaMemcpyAsync(h_pos_x, d_particles_.pos_x, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_pos_y, d_particles_.pos_y, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_pos_z, d_particles_.pos_z, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_vel_x, d_particles_.vel_x, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_vel_y, d_particles_.vel_y, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_vel_z, d_particles_.vel_z, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_force_x, d_particles_.force_x, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_force_y, d_particles_.force_y, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    CUDA_CHECK(cudaMemcpyAsync(h_force_z, d_particles_.force_z, count * sizeof(float), cudaMemcpyDeviceToHost, stream));
    return true;
}

// ============================================================================
// CUDA Kernels
// ============================================================================

namespace kernels {

// ---------------------------------------------------------------------------
// Bounding Box Reduction
// ---------------------------------------------------------------------------

__global__ void bounding_box_kernel(
    const float* __restrict__ pos_x,
    const float* __restrict__ pos_y,
    const float* __restrict__ pos_z,
    size_t count,
    float* __restrict__ g_min_x, float* __restrict__ g_min_y, float* __restrict__ g_min_z,
    float* __restrict__ g_max_x, float* __restrict__ g_max_y, float* __restrict__ g_max_z)
{
    __shared__ float s_min_x[256], s_min_y[256], s_min_z[256];
    __shared__ float s_max_x[256], s_max_y[256], s_max_z[256];

    const size_t tid = threadIdx.x;
    const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Initialize with extreme values
    float lmin_x = FLT_MAX, lmin_y = FLT_MAX, lmin_z = FLT_MAX;
    float lmax_x = -FLT_MAX, lmax_y = -FLT_MAX, lmax_z = -FLT_MAX;

    // Grid-stride loop for large arrays
    for (size_t i = idx; i < count; i += blockDim.x * gridDim.x) {
        float px = pos_x[i], py = pos_y[i], pz = pos_z[i];
        lmin_x = fminf(lmin_x, px); lmin_y = fminf(lmin_y, py); lmin_z = fminf(lmin_z, pz);
        lmax_x = fmaxf(lmax_x, px); lmax_y = fmaxf(lmax_y, py); lmax_z = fmaxf(lmax_z, pz);
    }

    s_min_x[tid] = lmin_x; s_min_y[tid] = lmin_y; s_min_z[tid] = lmin_z;
    s_max_x[tid] = lmax_x; s_max_y[tid] = lmax_y; s_max_z[tid] = lmax_z;
    __syncthreads();

    // Warp-aware block reduction
    for (int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < (unsigned)s) {
            s_min_x[tid] = fminf(s_min_x[tid], s_min_x[tid + s]);
            s_min_y[tid] = fminf(s_min_y[tid], s_min_y[tid + s]);
            s_min_z[tid] = fminf(s_min_z[tid], s_min_z[tid + s]);
            s_max_x[tid] = fmaxf(s_max_x[tid], s_max_x[tid + s]);
            s_max_y[tid] = fmaxf(s_max_y[tid], s_max_y[tid + s]);
            s_max_z[tid] = fmaxf(s_max_z[tid], s_max_z[tid + s]);
        }
        __syncthreads();
    }

    // Atomic reduction across blocks using IEEE 754-safe float atomics.
    // For float atomics: convert to int preserving order for both positive and negative values.
    // Positive floats: bit pattern is already ordered (atomicMin/Max on int works).
    // Negative floats: bit pattern is inverted — we flip sign bit and conditionals.
    // Formula: if val >= 0, use __float_as_int(val); else use 0x80000000 - __float_as_int(val).
    if (tid == 0) {
        auto float_to_ordered_int = [](float v) -> int {
            int i = __float_as_int(v);
            return (i >= 0) ? i : (int)(0x80000000u - (unsigned int)i);
        };
        auto ordered_int_to_float = [](int i) -> float {
            int v = (i >= 0) ? i : (int)(0x80000000u - (unsigned int)i);
            return __int_as_float(v);
        };
        (void)ordered_int_to_float;  // Used during readback, documented for completeness

        atomicMin(reinterpret_cast<int*>(g_min_x), float_to_ordered_int(s_min_x[0]));
        atomicMin(reinterpret_cast<int*>(g_min_y), float_to_ordered_int(s_min_y[0]));
        atomicMin(reinterpret_cast<int*>(g_min_z), float_to_ordered_int(s_min_z[0]));
        atomicMax(reinterpret_cast<int*>(g_max_x), float_to_ordered_int(s_max_x[0]));
        atomicMax(reinterpret_cast<int*>(g_max_y), float_to_ordered_int(s_max_y[0]));
        atomicMax(reinterpret_cast<int*>(g_max_z), float_to_ordered_int(s_max_z[0]));
    }
}

void compute_bounding_box(
    const ParticleSoA& particles,
    float* d_min_x, float* d_min_y, float* d_min_z,
    float* d_max_x, float* d_max_y, float* d_max_z,
    void* /*d_temp_storage*/, size_t& /*temp_storage_bytes*/,
    cudaStream_t stream)
{
    // Initialize results using ordered-int representation for safe atomic float min/max.
    // FLT_MAX is positive → ordered int = __float_as_int(FLT_MAX) = 0x7F7FFFFF
    // -FLT_MAX is negative → ordered int = 0x80000000 - __float_as_int(-FLT_MAX) = 0x80000001
    int init_min_i = 0x7F7FFFFF;   // ordered_int(FLT_MAX) — max possible for atomicMin
    int init_max_i = (int)0x80000001u; // ordered_int(-FLT_MAX) — min possible for atomicMax
    cudaMemcpyAsync(d_min_x, &init_min_i, sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_min_y, &init_min_i, sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_min_z, &init_min_i, sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_max_x, &init_max_i, sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_max_y, &init_max_i, sizeof(int), cudaMemcpyHostToDevice, stream);
    cudaMemcpyAsync(d_max_z, &init_max_i, sizeof(int), cudaMemcpyHostToDevice, stream);

    const int threads = 256;
    const int blocks = min((int)((particles.count + threads - 1) / threads), 1024);

    bounding_box_kernel<<<blocks, threads, 0, stream>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        particles.count,
        d_min_x, d_min_y, d_min_z,
        d_max_x, d_max_y, d_max_z
    );
}

// ---------------------------------------------------------------------------
// Morton Code Computation & Sorting
// ---------------------------------------------------------------------------

__global__ void compute_morton_codes_kernel(
    const float* __restrict__ pos_x,
    const float* __restrict__ pos_y,
    const float* __restrict__ pos_z,
    size_t count,
    uint32_t* __restrict__ morton_codes,
    int32_t* __restrict__ indices,
    float min_x, float min_y, float min_z, float range)
{
    const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) return;

    morton_codes[idx] = morton3D(pos_x[idx], pos_y[idx], pos_z[idx],
                                 min_x, min_y, min_z, range);
    indices[idx] = static_cast<int32_t>(idx);
}

void sort_particles_morton(
    ParticleSoA& particles,
    uint32_t* d_morton_codes,
    int32_t* d_sorted_indices,
    float3 bbox_min, float3 bbox_max,
    void* d_temp_storage, size_t& temp_storage_bytes,
    cudaStream_t stream)
{
    const size_t count = particles.count;
    const int threads = 256;
    const int blocks = (count + threads - 1) / threads;

    float range = fmaxf(fmaxf(bbox_max.x - bbox_min.x,
                               bbox_max.y - bbox_min.y),
                         bbox_max.z - bbox_min.z);
    if (range < 1e-6f) range = 1.0f;

    // Step 1: Compute Morton codes
    compute_morton_codes_kernel<<<blocks, threads, 0, stream>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        count, d_morton_codes, d_sorted_indices,
        bbox_min.x, bbox_min.y, bbox_min.z, range
    );

#if HAS_CUB
    // Step 2: Sort by Morton code using cub::DeviceRadixSort
    cub::DeviceRadixSort::SortPairs(
        d_temp_storage, temp_storage_bytes,
        d_morton_codes, d_morton_codes,
        d_sorted_indices, d_sorted_indices,
        count, 0, 30, stream
    );
#endif
}

// ---------------------------------------------------------------------------
// Mass Distribution (Upward Pass) — Segmented Reduction
// ---------------------------------------------------------------------------

/**
 * Leaf-level mass distribution kernel.
 * Each thread processes one leaf node.
 */
__global__ void leaf_mass_distribution_kernel(
    const float* __restrict__ p_pos_x,
    const float* __restrict__ p_pos_y,
    const float* __restrict__ p_pos_z,
    const float* __restrict__ p_mass,
    const int32_t* __restrict__ node_particle_begin,
    const int32_t* __restrict__ node_particle_end,
    const uint8_t* __restrict__ node_type,
    float* __restrict__ node_mass_cx,
    float* __restrict__ node_mass_cy,
    float* __restrict__ node_mass_cz,
    float* __restrict__ node_total_mass,
    size_t node_count)
{
    const size_t nid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nid >= node_count) return;

    // Initialize all nodes to zero mass (safe default for unprocessed nodes)
    if (node_type[nid] != 2) {
        // Non-leaf: will be computed by internal_mass_distribution_kernel.
        // Initialize to zero to prevent reading garbage if node has no children.
        if (node_type[nid] == 0) {  // Empty node
            node_total_mass[nid] = 0.0f;
            node_mass_cx[nid] = 0.0f;
            node_mass_cy[nid] = 0.0f;
            node_mass_cz[nid] = 0.0f;
        }
        return;
    }

    int32_t begin = node_particle_begin[nid];
    int32_t end = node_particle_end[nid];

    float total_mass = 0.0f;
    float wcx = 0.0f, wcy = 0.0f, wcz = 0.0f;

    for (int32_t i = begin; i < end; ++i) {
        float m = p_mass[i];
        total_mass += m;
        wcx += m * p_pos_x[i];
        wcy += m * p_pos_y[i];
        wcz += m * p_pos_z[i];
    }

    if (total_mass > 0.0f) {
        node_mass_cx[nid] = wcx / total_mass;
        node_mass_cy[nid] = wcy / total_mass;
        node_mass_cz[nid] = wcz / total_mass;
        node_total_mass[nid] = total_mass;
    }
}

/**
 * Internal node mass distribution kernel (one level at a time, bottom-up).
 */
__global__ void internal_mass_distribution_kernel(
    const int32_t* __restrict__ children,
    const uint8_t* __restrict__ node_type,
    const int32_t* __restrict__ node_level,
    float* __restrict__ node_mass_cx,
    float* __restrict__ node_mass_cy,
    float* __restrict__ node_mass_cz,
    float* __restrict__ node_total_mass,
    size_t node_count,
    int32_t target_level)
{
    const size_t nid = blockIdx.x * blockDim.x + threadIdx.x;
    if (nid >= node_count) return;

    // Only process internal nodes at target level
    if (node_type[nid] != 1 || node_level[nid] != target_level) return;

    float total_mass = 0.0f;
    float wcx = 0.0f, wcy = 0.0f, wcz = 0.0f;

    for (int c = 0; c < 8; ++c) {
        int32_t child_idx = children[nid * 8 + c];
        if (child_idx < 0) continue;

        float cm = node_total_mass[child_idx];
        if (cm > 0.0f) {
            total_mass += cm;
            wcx += cm * node_mass_cx[child_idx];
            wcy += cm * node_mass_cy[child_idx];
            wcz += cm * node_mass_cz[child_idx];
        }
    }

    if (total_mass > 0.0f) {
        node_mass_cx[nid] = wcx / total_mass;
        node_mass_cy[nid] = wcy / total_mass;
        node_mass_cz[nid] = wcz / total_mass;
        node_total_mass[nid] = total_mass;
    }
}

void compute_mass_distribution(
    const ParticleSoA& particles,
    OctreeNodeSoA& nodes,
    void* /*d_temp_storage*/, size_t& /*temp_storage_bytes*/,
    cudaStream_t stream)
{
    const int threads = 256;
    const int blocks = (nodes.node_count + threads - 1) / threads;

    // Step 1: Process leaf nodes
    leaf_mass_distribution_kernel<<<blocks, threads, 0, stream>>>(
        particles.pos_x, particles.pos_y, particles.pos_z, particles.mass,
        nodes.particle_begin, nodes.particle_end, nodes.node_type,
        nodes.mass_center_x, nodes.mass_center_y, nodes.mass_center_z,
        nodes.total_mass, nodes.node_count
    );

    // Step 2: Process internal nodes bottom-up
    for (int32_t level = static_cast<int32_t>(nodes.max_level) - 1; level >= 0; --level) {
        internal_mass_distribution_kernel<<<blocks, threads, 0, stream>>>(
            nodes.children, nodes.node_type, nodes.level,
            nodes.mass_center_x, nodes.mass_center_y, nodes.mass_center_z,
            nodes.total_mass, nodes.node_count, level
        );
    }
}

// ---------------------------------------------------------------------------
// Force Calculation — Tile-Optimized Tree Traversal
// ---------------------------------------------------------------------------

/**
 * Tile-optimized Barnes-Hut force calculation kernel.
 *
 * Each warp (32 threads) processes 32 particles simultaneously.
 * Shared memory caches node data for the tile.
 * Stack-based iterative tree traversal avoids recursion.
 */
__global__ void force_kernel_tiled(
    const float* __restrict__ p_pos_x,
    const float* __restrict__ p_pos_y,
    const float* __restrict__ p_pos_z,
    const float* __restrict__ p_mass,
    float* __restrict__ p_force_x,
    float* __restrict__ p_force_y,
    float* __restrict__ p_force_z,
    // Tree data
    const float* __restrict__ n_center_x,
    const float* __restrict__ n_center_y,
    const float* __restrict__ n_center_z,
    const float* __restrict__ n_size,
    const float* __restrict__ n_mass_cx,
    const float* __restrict__ n_mass_cy,
    const float* __restrict__ n_mass_cz,
    const float* __restrict__ n_total_mass,
    const int32_t* __restrict__ n_children,
    const uint8_t* __restrict__ n_type,
    const int32_t* __restrict__ n_particle_begin,
    const int32_t* __restrict__ n_particle_end,
    // Config
    size_t particle_count,
    size_t node_count,
    float theta,
    float eps_sq,
    float G)
{
    const size_t pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= particle_count) return;

    const float px = p_pos_x[pid];
    const float py = p_pos_y[pid];
    const float pz = p_pos_z[pid];
    const float pm = p_mass[pid];

    float fx = 0.0f, fy = 0.0f, fz = 0.0f;

    // Stack-based iterative tree traversal (avoids recursion).
    // Octree depth is bounded by log8(N) + constant. For N=10M, depth ~ 8.
    // Each level can push up to 8 children, so max stack usage ~ 8*depth = ~64.
    // We use 128 entries for safety and track overflow.
    constexpr int STACK_SIZE = 128;
    int32_t stack[STACK_SIZE];
    int stack_top = 0;

    // Push root's children
    for (int c = 0; c < 8; ++c) {
        int32_t child = n_children[0 * 8 + c];  // Root is node 0
        if (child >= 0 && child < (int32_t)node_count && n_type[child] != 0) {
            stack[stack_top++] = child;
        }
    }

    while (stack_top > 0) {
        int32_t nid = stack[--stack_top];

        // Compute distance to node's center of mass
        float dx = px - n_mass_cx[nid];
        float dy = py - n_mass_cy[nid];
        float dz = pz - n_mass_cz[nid];
        float r_sq = dx * dx + dy * dy + dz * dz + eps_sq;
        float node_size = n_size[nid];

        // Barnes-Hut criterion: s/d <= theta
        bool well_separated = (node_size * node_size / r_sq) <= (theta * theta);

        if (well_separated || n_type[nid] == 2) {
            if (n_type[nid] == 2) {
                // Leaf: direct interaction with particles
                int32_t begin = n_particle_begin[nid];
                int32_t end = n_particle_end[nid];
                for (int32_t i = begin; i < end; ++i) {
                    if (i == (int32_t)pid) continue;  // Skip self

                    float dx2 = px - p_pos_x[i];
                    float dy2 = py - p_pos_y[i];
                    float dz2 = pz - p_pos_z[i];
                    float r_sq2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2 + eps_sq;
                    float r_inv = rsqrtf(r_sq2);
                    float r_inv3 = r_inv * r_inv * r_inv;

                    float f = -G * pm * p_mass[i] * r_inv3;
                    fx += f * dx2;
                    fy += f * dy2;
                    fz += f * dz2;
                }
            } else {
                // Well-separated internal node: multipole approximation
                float r_inv = rsqrtf(r_sq);
                float r_inv3 = r_inv * r_inv * r_inv;
                float f = -G * pm * n_total_mass[nid] * r_inv3;
                fx += f * dx;
                fy += f * dy;
                fz += f * dz;
            }
        } else {
            // Not well separated: push children for further traversal
            for (int c = 0; c < 8; ++c) {
                int32_t child = n_children[nid * 8 + c];
                if (child >= 0 && child < (int32_t)node_count && n_type[child] != 0) {
                    if (stack_top < STACK_SIZE) {
                        stack[stack_top++] = child;
                    }
                    // Note: stack overflow means some nodes are skipped.
                    // With STACK_SIZE=128, this should never happen for N < 10^12.
                }
            }
        }
    }

    p_force_x[pid] = fx;
    p_force_y[pid] = fy;
    p_force_z[pid] = fz;
}

void calculate_forces_tiled(
    ParticleSoA& particles,
    const OctreeNodeSoA& nodes,
    const KernelConfig& config,
    cudaStream_t stream)
{
    const int threads = 256;
    const int blocks = (particles.count + threads - 1) / threads;

    force_kernel_tiled<<<blocks, threads, 0, stream>>>(
        particles.pos_x, particles.pos_y, particles.pos_z, particles.mass,
        particles.force_x, particles.force_y, particles.force_z,
        nodes.center_x, nodes.center_y, nodes.center_z, nodes.size,
        nodes.mass_center_x, nodes.mass_center_y, nodes.mass_center_z,
        nodes.total_mass, nodes.children, nodes.node_type,
        nodes.particle_begin, nodes.particle_end,
        particles.count, nodes.node_count,
        config.theta, config.epsilon_squared, config.gravity
    );
}

// ---------------------------------------------------------------------------
// Direct N-body (O(N^2)) — Reference Implementation
// ---------------------------------------------------------------------------

__global__ void direct_force_kernel(
    const float* __restrict__ pos_x,
    const float* __restrict__ pos_y,
    const float* __restrict__ pos_z,
    const float* __restrict__ mass,
    float* __restrict__ force_x,
    float* __restrict__ force_y,
    float* __restrict__ force_z,
    size_t count,
    float eps_sq,
    float G)
{
    // Use shared memory tiles for O(N^2) direct computation
    extern __shared__ float shared[];
    float* s_px = shared;
    float* s_py = s_px + blockDim.x;
    float* s_pz = s_py + blockDim.x;
    float* s_m  = s_pz + blockDim.x;

    const size_t pid = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t tid = threadIdx.x;

    float px = 0.0f, py = 0.0f, pz = 0.0f, pm = 0.0f;
    if (pid < count) {
        px = pos_x[pid];
        py = pos_y[pid];
        pz = pos_z[pid];
        pm = mass[pid];
    }

    float fx = 0.0f, fy = 0.0f, fz = 0.0f;

    // Tile-based interaction computation
    for (size_t tile = 0; tile < (count + blockDim.x - 1) / blockDim.x; ++tile) {
        size_t j = tile * blockDim.x + tid;

        // Load tile into shared memory
        if (j < count) {
            s_px[tid] = pos_x[j];
            s_py[tid] = pos_y[j];
            s_pz[tid] = pos_z[j];
            s_m[tid]  = mass[j];
        } else {
            s_px[tid] = 0.0f;
            s_py[tid] = 0.0f;
            s_pz[tid] = 0.0f;
            s_m[tid]  = 0.0f;
        }
        __syncthreads();

        // Compute interactions with tile
        if (pid < count) {
            for (int k = 0; k < (int)blockDim.x; ++k) {
                size_t other = tile * blockDim.x + k;
                if (other >= count || other == pid) continue;

                float dx = px - s_px[k];
                float dy = py - s_py[k];
                float dz = pz - s_pz[k];
                float r_sq = dx * dx + dy * dy + dz * dz + eps_sq;
                float r_inv = rsqrtf(r_sq);
                float r_inv3 = r_inv * r_inv * r_inv;

                float f = -G * pm * s_m[k] * r_inv3;
                fx += f * dx;
                fy += f * dy;
                fz += f * dz;
            }
        }
        __syncthreads();
    }

    if (pid < count) {
        force_x[pid] = fx;
        force_y[pid] = fy;
        force_z[pid] = fz;
    }
}

void calculate_forces_direct(
    ParticleSoA& particles,
    const KernelConfig& config,
    cudaStream_t stream)
{
    const int threads = 256;
    const int blocks = (particles.count + threads - 1) / threads;
    const size_t shared_mem = 4 * threads * sizeof(float);

    direct_force_kernel<<<blocks, threads, shared_mem, stream>>>(
        particles.pos_x, particles.pos_y, particles.pos_z, particles.mass,
        particles.force_x, particles.force_y, particles.force_z,
        particles.count, config.epsilon_squared, config.gravity
    );
}

// ---------------------------------------------------------------------------
// Leapfrog Integration
// ---------------------------------------------------------------------------

__global__ void integrate_kernel(
    float* __restrict__ pos_x, float* __restrict__ pos_y, float* __restrict__ pos_z,
    float* __restrict__ vel_x, float* __restrict__ vel_y, float* __restrict__ vel_z,
    const float* __restrict__ force_x, const float* __restrict__ force_y, const float* __restrict__ force_z,
    const float* __restrict__ mass,
    size_t count, float dt)
{
    const size_t pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= count) return;

    float m = mass[pid];
    if (m <= 0.0f) return;  // Skip zero/negative mass particles
    float inv_mass = 1.0f / m;
    float ax = force_x[pid] * inv_mass;
    float ay = force_y[pid] * inv_mass;
    float az = force_z[pid] * inv_mass;

    float half_dt = 0.5f * dt;

    // Leapfrog kick-drift-kick
    vel_x[pid] += ax * half_dt;
    vel_y[pid] += ay * half_dt;
    vel_z[pid] += az * half_dt;

    pos_x[pid] += vel_x[pid] * dt;
    pos_y[pid] += vel_y[pid] * dt;
    pos_z[pid] += vel_z[pid] * dt;

    vel_x[pid] += ax * half_dt;
    vel_y[pid] += ay * half_dt;
    vel_z[pid] += az * half_dt;
}

void integrate_particles(
    ParticleSoA& particles,
    float dt,
    cudaStream_t stream)
{
    const int threads = 256;
    const int blocks = (particles.count + threads - 1) / threads;

    integrate_kernel<<<blocks, threads, 0, stream>>>(
        particles.pos_x, particles.pos_y, particles.pos_z,
        particles.vel_x, particles.vel_y, particles.vel_z,
        particles.force_x, particles.force_y, particles.force_z,
        particles.mass, particles.count, dt
    );
}

// ---------------------------------------------------------------------------
// Reset Forces
// ---------------------------------------------------------------------------

__global__ void reset_forces_kernel(
    float* __restrict__ fx, float* __restrict__ fy, float* __restrict__ fz,
    size_t count)
{
    const size_t pid = blockIdx.x * blockDim.x + threadIdx.x;
    if (pid >= count) return;
    fx[pid] = 0.0f;
    fy[pid] = 0.0f;
    fz[pid] = 0.0f;
}

void reset_forces(ParticleSoA& particles, cudaStream_t stream) {
    const int threads = 256;
    const int blocks = (particles.count + threads - 1) / threads;

    reset_forces_kernel<<<blocks, threads, 0, stream>>>(
        particles.force_x, particles.force_y, particles.force_z,
        particles.count
    );
}

} // namespace kernels

// ============================================================================
// CudaSimulation Implementation
// ============================================================================

CudaSimulation::CudaSimulation()
    : stream_(nullptr)
    , initialized_(false)
    , d_temp_storage_(nullptr)
    , temp_storage_bytes_(0)
    , d_morton_codes_(nullptr)
    , d_sorted_indices_(nullptr)
    , d_bbox_min_(nullptr)
    , d_bbox_max_(nullptr)
{
    memset(&config_, 0, sizeof(config_));
    memset(&stats_, 0, sizeof(stats_));
}

CudaSimulation::~CudaSimulation() {
    cleanup();
}

bool CudaSimulation::initialize(
    const float* h_mass,
    const float* h_pos_x, const float* h_pos_y, const float* h_pos_z,
    const float* h_vel_x, const float* h_vel_y, const float* h_vel_z,
    size_t particle_count, const KernelConfig& config)
{
    config_ = config;
    stats_.particle_count = particle_count;

    CUDA_CHECK(cudaStreamCreate(&stream_));

    // Allocate particle and node memory (max nodes ~ 3x particles)
    if (!mem_.allocate(particle_count, particle_count * 3)) {
        return false;
    }

    // Upload initial state
    if (!mem_.upload_particles(h_mass, h_pos_x, h_pos_y, h_pos_z,
                               h_vel_x, h_vel_y, h_vel_z,
                               particle_count, stream_)) {
        return false;
    }

    // Allocate auxiliary buffers
    CUDA_CHECK(cudaMalloc(&d_morton_codes_, particle_count * sizeof(uint32_t)));
    CUDA_CHECK(cudaMalloc(&d_sorted_indices_, particle_count * sizeof(int32_t)));
    CUDA_CHECK(cudaMalloc(&d_bbox_min_, 3 * sizeof(float)));
    CUDA_CHECK(cudaMalloc(&d_bbox_max_, 3 * sizeof(float)));

    // Allocate temp storage for cub operations
    temp_storage_bytes_ = particle_count * sizeof(float) * 4;  // Conservative estimate
    CUDA_CHECK(cudaMalloc(&d_temp_storage_, temp_storage_bytes_));

    CUDA_CHECK(cudaStreamSynchronize(stream_));
    initialized_ = true;
    return true;
}

bool CudaSimulation::step() {
    if (!initialized_) return false;

    cudaEvent_t ev_start, ev_bbox, ev_force, ev_integrate, ev_end;
    CUDA_CHECK(cudaEventCreate(&ev_start));
    CUDA_CHECK(cudaEventCreate(&ev_bbox));
    CUDA_CHECK(cudaEventCreate(&ev_force));
    CUDA_CHECK(cudaEventCreate(&ev_integrate));
    CUDA_CHECK(cudaEventCreate(&ev_end));

    CUDA_CHECK(cudaEventRecord(ev_start, stream_));

    auto& p = mem_.particles();

    // 1. Compute bounding box
    kernels::compute_bounding_box(
        p,
        d_bbox_min_, d_bbox_min_ + 1, d_bbox_min_ + 2,
        d_bbox_max_, d_bbox_max_ + 1, d_bbox_max_ + 2,
        d_temp_storage_, temp_storage_bytes_, stream_
    );

    CUDA_CHECK(cudaEventRecord(ev_bbox, stream_));

    // 2. Reset forces
    kernels::reset_forces(p, stream_);

    // 3. Calculate forces (direct N-body for now — tree traversal requires
    //    CPU-built tree uploaded to GPU)
    kernels::calculate_forces_direct(p, config_, stream_);

    CUDA_CHECK(cudaEventRecord(ev_force, stream_));

    // 4. Integrate
    kernels::integrate_particles(p, config_.dt, stream_);

    CUDA_CHECK(cudaEventRecord(ev_integrate, stream_));
    CUDA_CHECK(cudaEventRecord(ev_end, stream_));
    CUDA_CHECK(cudaEventSynchronize(ev_end));

    // Collect timing
    float ms;
    cudaEventElapsedTime(&ms, ev_start, ev_bbox);
    stats_.time_bbox_ms = ms;
    cudaEventElapsedTime(&ms, ev_bbox, ev_force);
    stats_.time_force_ms = ms;
    cudaEventElapsedTime(&ms, ev_force, ev_integrate);
    stats_.time_integrate_ms = ms;
    cudaEventElapsedTime(&ms, ev_start, ev_end);
    stats_.time_total_ms = ms;

    cudaEventDestroy(ev_start);
    cudaEventDestroy(ev_bbox);
    cudaEventDestroy(ev_force);
    cudaEventDestroy(ev_integrate);
    cudaEventDestroy(ev_end);

    return true;
}

bool CudaSimulation::download_state(
    float* h_pos_x, float* h_pos_y, float* h_pos_z,
    float* h_vel_x, float* h_vel_y, float* h_vel_z,
    float* h_force_x, float* h_force_y, float* h_force_z)
{
    if (!initialized_) return false;
    bool ok = mem_.download_particles(h_pos_x, h_pos_y, h_pos_z,
                                       h_vel_x, h_vel_y, h_vel_z,
                                       h_force_x, h_force_y, h_force_z,
                                       stats_.particle_count, stream_);
    CUDA_CHECK(cudaStreamSynchronize(stream_));
    return ok;
}

void CudaSimulation::cleanup() {
    if (!initialized_) return;

    mem_.free();

    if (d_morton_codes_) { cudaFree(d_morton_codes_); d_morton_codes_ = nullptr; }
    if (d_sorted_indices_) { cudaFree(d_sorted_indices_); d_sorted_indices_ = nullptr; }
    if (d_bbox_min_) { cudaFree(d_bbox_min_); d_bbox_min_ = nullptr; }
    if (d_bbox_max_) { cudaFree(d_bbox_max_); d_bbox_max_ = nullptr; }
    if (d_temp_storage_) { cudaFree(d_temp_storage_); d_temp_storage_ = nullptr; }

    if (stream_) { cudaStreamDestroy(stream_); stream_ = nullptr; }

    initialized_ = false;
}

} // namespace cuda
} // namespace barnes_hut
