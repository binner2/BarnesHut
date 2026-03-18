#pragma once

/**
 * Barnes-Hut CUDA 13.2 Enhanced Kernels
 *
 * GPU-accelerated core algorithm using CUDA 13.2 features:
 * - Structure-of-Arrays (SoA) data layout for coalesced memory access
 * - cub::DeviceReduce for bounding box
 * - cub::DeviceSegmentedReduce for mass distribution
 * - Tile-optimized force calculation kernel
 * - CUDA Graphs for simulation loop
 * - Modern CCCL 3.2 API usage
 */

#include <cuda_runtime.h>
#include <cstddef>
#include <cstdint>

namespace barnes_hut {
namespace cuda {

// ============================================================================
// SoA (Structure-of-Arrays) Data Layout for GPU
// ============================================================================

/**
 * GPU-friendly particle data in SoA format.
 * Each field is a separate contiguous array for coalesced memory access.
 */
struct ParticleSoA {
    // Positions (read/write)
    float* pos_x;
    float* pos_y;
    float* pos_z;

    // Velocities (read/write)
    float* vel_x;
    float* vel_y;
    float* vel_z;

    // Forces (write during force calc, read during integration)
    float* force_x;
    float* force_y;
    float* force_z;

    // Mass (read-only after initialization)
    float* mass;

    // Particle count
    size_t count;
};

/**
 * GPU-friendly octree node in SoA format (linearized).
 * Nodes are stored level-by-level for bottom-up processing.
 */
struct OctreeNodeSoA {
    // Geometric properties
    float* center_x;
    float* center_y;
    float* center_z;
    float* size;

    // Mass distribution (computed in upward pass)
    float* mass_center_x;
    float* mass_center_y;
    float* mass_center_z;
    float* total_mass;

    // Tree structure (indices, not pointers)
    int32_t* children;         // [node_count * 8] — child indices (-1 = empty)
    int32_t* parent;           // Parent index
    int32_t* particle_begin;   // First particle index in this node
    int32_t* particle_end;     // One-past-last particle index
    uint8_t* node_type;        // 0=empty, 1=internal, 2=leaf
    int32_t* level;            // Depth level

    // Metadata
    size_t node_count;
    size_t max_level;
};

// ============================================================================
// Morton Code Utilities
// ============================================================================

/**
 * Compute 30-bit Morton code for 3D position.
 * Used for spatial sorting to improve memory locality.
 */
__device__ __host__ inline uint32_t morton3D(float x, float y, float z,
                                              float min_x, float min_y, float min_z,
                                              float range) {
    // Normalize to [0, 1023], clamp BEFORE cast to prevent undefined behavior
    float nx = fminf(fmaxf(((x - min_x) / range), 0.0f), 1.0f) * 1023.0f;
    float ny = fminf(fmaxf(((y - min_y) / range), 0.0f), 1.0f) * 1023.0f;
    float nz = fminf(fmaxf(((z - min_z) / range), 0.0f), 1.0f) * 1023.0f;
    uint32_t ix = static_cast<uint32_t>(nx);
    uint32_t iy = static_cast<uint32_t>(ny);
    uint32_t iz = static_cast<uint32_t>(nz);

    // Spread bits: insert two 0-bits after each bit
    auto spread = [](uint32_t v) -> uint32_t {
        v = (v | (v << 16)) & 0x030000FF;
        v = (v | (v <<  8)) & 0x0300F00F;
        v = (v | (v <<  4)) & 0x030C30C3;
        v = (v | (v <<  2)) & 0x09249249;
        return v;
    };

    return spread(ix) | (spread(iy) << 1) | (spread(iz) << 2);
}

// ============================================================================
// Kernel Configuration
// ============================================================================

struct KernelConfig {
    float theta;                // Barnes-Hut opening angle
    float epsilon_squared;      // Gravitational softening
    float gravity;              // Gravitational constant
    float dt;                   // Time step
    int max_particles_per_leaf; // Leaf capacity
};

// ============================================================================
// GPU Memory Manager
// ============================================================================

/**
 * Manages all GPU memory for the simulation.
 * RAII-based allocation and deallocation.
 */
class GpuMemoryManager {
public:
    GpuMemoryManager();
    ~GpuMemoryManager();

    // Non-copyable
    GpuMemoryManager(const GpuMemoryManager&) = delete;
    GpuMemoryManager& operator=(const GpuMemoryManager&) = delete;

    /**
     * Allocate GPU memory for particles and tree.
     * @param particle_count Number of particles
     * @param max_node_count Maximum number of tree nodes (typically 3*N)
     * @return true if allocation succeeded
     */
    bool allocate(size_t particle_count, size_t max_node_count);

    /**
     * Free all GPU memory.
     */
    void free();

    /**
     * Upload particle data from host AoS to device SoA.
     * Performs AoS→SoA conversion on the fly.
     */
    bool upload_particles(const float* h_mass,
                          const float* h_pos_x, const float* h_pos_y, const float* h_pos_z,
                          const float* h_vel_x, const float* h_vel_y, const float* h_vel_z,
                          size_t count, cudaStream_t stream);

    /**
     * Download particle results from device SoA to host.
     */
    bool download_particles(float* h_pos_x, float* h_pos_y, float* h_pos_z,
                            float* h_vel_x, float* h_vel_y, float* h_vel_z,
                            float* h_force_x, float* h_force_y, float* h_force_z,
                            size_t count, cudaStream_t stream);

    // Accessors
    ParticleSoA& particles() { return d_particles_; }
    OctreeNodeSoA& nodes() { return d_nodes_; }
    const ParticleSoA& particles() const { return d_particles_; }
    const OctreeNodeSoA& nodes() const { return d_nodes_; }

private:
    ParticleSoA d_particles_;
    OctreeNodeSoA d_nodes_;
    bool allocated_;
};

// ============================================================================
// Kernel Launch Functions
// ============================================================================

namespace kernels {

/**
 * Compute bounding box using parallel reduction.
 * Uses cub::DeviceReduce::Min/Max internally.
 *
 * @param particles SoA particle data
 * @param[out] min_x/y/z Bounding box minimum
 * @param[out] max_x/y/z Bounding box maximum
 * @param stream CUDA stream
 */
void compute_bounding_box(
    const ParticleSoA& particles,
    float* d_min_x, float* d_min_y, float* d_min_z,
    float* d_max_x, float* d_max_y, float* d_max_z,
    void* d_temp_storage, size_t& temp_storage_bytes,
    cudaStream_t stream
);

/**
 * Sort particles by Morton code for spatial locality.
 *
 * @param particles SoA particle data (reordered in place)
 * @param d_morton_codes Device array for Morton codes
 * @param d_sorted_indices Device array for sorted indices
 * @param bbox_min/max Bounding box
 * @param stream CUDA stream
 */
void sort_particles_morton(
    ParticleSoA& particles,
    uint32_t* d_morton_codes,
    int32_t* d_sorted_indices,
    float3 bbox_min, float3 bbox_max,
    void* d_temp_storage, size_t& temp_storage_bytes,
    cudaStream_t stream
);

/**
 * Compute mass distribution (upward pass) using segmented reduction.
 * Uses cub::DeviceSegmentedReduce for leaf-level computation
 * and level-by-level kernel for internal nodes.
 *
 * @param particles SoA particle data
 * @param nodes SoA node data
 * @param stream CUDA stream
 */
void compute_mass_distribution(
    const ParticleSoA& particles,
    OctreeNodeSoA& nodes,
    void* d_temp_storage, size_t& temp_storage_bytes,
    cudaStream_t stream
);

/**
 * Calculate gravitational forces using tile-optimized traversal.
 * Each tile = 32 particles (warp-aligned).
 * Shared memory used for node data caching.
 *
 * @param particles SoA particle data
 * @param nodes SoA node data (read-only)
 * @param config Simulation parameters
 * @param stream CUDA stream
 */
void calculate_forces_tiled(
    ParticleSoA& particles,
    const OctreeNodeSoA& nodes,
    const KernelConfig& config,
    cudaStream_t stream
);

/**
 * Direct N-body force calculation (brute force, O(N^2)).
 * Used as reference for correctness validation.
 *
 * @param particles SoA particle data
 * @param config Simulation parameters
 * @param stream CUDA stream
 */
void calculate_forces_direct(
    ParticleSoA& particles,
    const KernelConfig& config,
    cudaStream_t stream
);

/**
 * Leapfrog integration kernel.
 * Trivially parallel — each particle independent.
 *
 * @param particles SoA particle data
 * @param dt Time step
 * @param stream CUDA stream
 */
void integrate_particles(
    ParticleSoA& particles,
    float dt,
    cudaStream_t stream
);

/**
 * Reset force accumulators to zero.
 *
 * @param particles SoA particle data
 * @param stream CUDA stream
 */
void reset_forces(
    ParticleSoA& particles,
    cudaStream_t stream
);

} // namespace kernels

// ============================================================================
// Simulation Runner
// ============================================================================

/**
 * High-level simulation runner using CUDA 13.2 features.
 * Manages the full simulation pipeline on GPU.
 */
class CudaSimulation {
public:
    struct Statistics {
        double time_bbox_ms;
        double time_sort_ms;
        double time_tree_ms;
        double time_mass_ms;
        double time_force_ms;
        double time_integrate_ms;
        double time_total_ms;
        size_t particle_count;
        size_t node_count;
    };

    CudaSimulation();
    ~CudaSimulation();

    /**
     * Initialize with particle data.
     * Allocates GPU memory and uploads initial state.
     */
    bool initialize(const float* h_mass,
                    const float* h_pos_x, const float* h_pos_y, const float* h_pos_z,
                    const float* h_vel_x, const float* h_vel_y, const float* h_vel_z,
                    size_t particle_count, const KernelConfig& config);

    /**
     * Run one simulation step on GPU.
     */
    bool step();

    /**
     * Download current state to host.
     */
    bool download_state(float* h_pos_x, float* h_pos_y, float* h_pos_z,
                        float* h_vel_x, float* h_vel_y, float* h_vel_z,
                        float* h_force_x, float* h_force_y, float* h_force_z);

    const Statistics& get_statistics() const { return stats_; }

    void cleanup();

private:
    GpuMemoryManager mem_;
    KernelConfig config_;
    Statistics stats_;
    cudaStream_t stream_;
    bool initialized_;

    // Temporary storage for cub algorithms
    void* d_temp_storage_;
    size_t temp_storage_bytes_;

    // Morton sort buffers
    uint32_t* d_morton_codes_;
    int32_t* d_sorted_indices_;

    // Bounding box results
    float* d_bbox_min_;  // [3]
    float* d_bbox_max_;  // [3]
};

} // namespace cuda
} // namespace barnes_hut
