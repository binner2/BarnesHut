#pragma once

/**
 * GPU Barnes-Hut N-body simulation (CUDA 13.4 toolchain, CCCL 3.4).
 *
 * Pipeline per force evaluation, entirely on the device:
 *   bounding box   cub::DeviceReduce over a cuda::transform_iterator
 *   Morton keys    63-bit (21 bits/axis)
 *   sort           cub::DeviceRadixSort::SortPairs (stable)
 *   gather         particle state permuted into key order
 *   hierarchy      Karras (2012) binary radix tree, one thread per node
 *   refit          bottom-up, cuda::atomic_ref visit counters (acq_rel)
 *   forces         per-particle stack traversal, monopole + softening
 * Time integration is kick-drift-kick leapfrog. The step can optionally be
 * captured once into a CUDA Graph and replayed (static sizes, preallocated
 * temporary storage), which removes per-kernel launch overhead.
 *
 * Particle data stays in Morton order on the device between steps; the
 * download functions return it in the caller's original order.
 */

#include "bh_core.cuh"

#include <cuda_runtime.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>

namespace barnes_hut {
namespace gpu {

class CudaError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct GpuConfig {
    float theta = 0.5f;               // opening parameter (size / distance)
    float epsilon_squared = 1e-10f;   // Plummer softening, squared
    float gravity = 1.0f;
    float dt = 1e-3f;
    bool use_graph = false;           // capture step() into a CUDA Graph
    bool collect_stats = false;       // count interactions (small overhead)
};

// Milliseconds measured with CUDA events for the most recent call.
struct GpuTiming {
    float tree_ms = 0.0f;       // bbox + keys + sort + gather + build + refit
    float force_ms = 0.0f;
    float integrate_ms = 0.0f;  // both kicks + drift
    float total_ms = 0.0f;
};

struct GpuForceStats {
    unsigned long long particle_particle = 0;
    unsigned long long particle_cell = 0;
    unsigned long long stack_overflows = 0;  // > 0 means degraded accuracy
};

struct CudaFree {
    void operator()(void* p) const noexcept { cudaFree(p); }
};

template <typename T>
using DeviceBuffer = std::unique_ptr<T[], CudaFree>;

class BarnesHutGpu {
public:
    // Uploads the initial state and computes the initial forces.
    BarnesHutGpu(const float* mass,
                 const float* pos_x, const float* pos_y, const float* pos_z,
                 const float* vel_x, const float* vel_y, const float* vel_z,
                 std::size_t count, const GpuConfig& config);
    ~BarnesHutGpu();

    BarnesHutGpu(const BarnesHutGpu&) = delete;
    BarnesHutGpu& operator=(const BarnesHutGpu&) = delete;

    // Barnes-Hut forces at the current positions (no time advance).
    void compute_forces();

    // O(N^2) direct-summation forces at the current positions (reference).
    void compute_forces_direct();

    // One kick-drift-kick step using Barnes-Hut forces.
    void step();

    // Wait for all queued work; returns timing of the last call.
    GpuTiming synchronize();

    // Device -> host, in the original particle order. Any pointer may be null.
    void download(float* pos_x, float* pos_y, float* pos_z,
                  float* vel_x, float* vel_y, float* vel_z,
                  float* force_x, float* force_y, float* force_z);

    [[nodiscard]] GpuForceStats force_stats();
    [[nodiscard]] std::size_t size() const noexcept { return n_; }

    static bool device_available() noexcept;

private:
    void enqueue_tree();
    void enqueue_forces();
    void enqueue_step(bool timed);
    void record(int slot);
    void capture_step_graph();

    std::size_t n_;
    GpuConfig config_;
    cudaStream_t stream_ = nullptr;

    // Particle state (Morton order) and gather scratch
    DeviceBuffer<float> px_, py_, pz_, vx_, vy_, vz_, m_, fx_, fy_, fz_;
    DeviceBuffer<float> tmp_px_, tmp_py_, tmp_pz_, tmp_vx_, tmp_vy_, tmp_vz_, tmp_m_;
    DeviceBuffer<int> id_, tmp_id_;  // sorted position -> original index

    // Sorting
    DeviceBuffer<std::uint64_t> keys_in_, keys_out_;
    DeviceBuffer<int> perm_in_, perm_out_;

    // Tree
    DeviceBuffer<int> left_, right_, parent_, visits_;
    DeviceBuffer<float4> cm_, lo_, hi_;

    // Reductions / diagnostics
    DeviceBuffer<unsigned char> temp_storage_;
    std::size_t temp_storage_bytes_ = 0;
    DeviceBuffer<Box> box_;
    DeviceBuffer<unsigned long long> counters_;     // [pp, pc, overflow]

    // Launch geometry (grid-stride kernels, from occupancy)
    int grid_elem_ = 1;
    int grid_force_ = 1;

    // Timing events: start, after kick+drift, after tree, after force, end.
    // In graph mode only start/end are recorded (total time).
    cudaEvent_t events_[5] = {};
    bool phase_timing_valid_ = false;

    cudaGraph_t graph_ = nullptr;
    cudaGraphExec_t graph_exec_ = nullptr;
};

} // namespace gpu
} // namespace barnes_hut
