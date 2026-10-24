/**
 * GPU Barnes-Hut N-body simulation — kernels and host driver.
 *
 * Kernels are grid-stride wrappers around the __host__ __device__ functions
 * in bh_core.cuh (verified on the host by tests/host_core_test.cu). Nothing
 * here requires a particular architecture beyond the toolkit's minimum
 * (sm_75 for CUDA 13.x): no clusters, no PDL, no TMA.
 */

#include "bh_gpu.cuh"

#include <cub/block/block_reduce.cuh>
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_reduce.cuh>
#include <cuda/iterator>

#include <algorithm>
#include <cmath>
#include <vector>

namespace barnes_hut {
namespace gpu {

namespace {

void check(cudaError_t err, const char* what, const char* file, int line) {
    if (err != cudaSuccess) {
        throw CudaError(std::string(what) + " failed at " + file + ":" + std::to_string(line) +
                        ": " + cudaGetErrorString(err));
    }
}

#define BH_CUDA_TRY(expr) check((expr), #expr, __FILE__, __LINE__)
#define BH_CUDA_LAUNCH_CHECK() check(cudaGetLastError(), "kernel launch", __FILE__, __LINE__)

template <typename T>
DeviceBuffer<T> device_alloc(std::size_t count) {
    T* ptr = nullptr;
    BH_CUDA_TRY(cudaMalloc(&ptr, std::max<std::size_t>(count, 1) * sizeof(T)));
    return DeviceBuffer<T>(ptr);
}

constexpr int ELEM_BLOCK = 256;
constexpr int FORCE_BLOCK = 128;
constexpr int TRAVERSAL_STACK = 64;  // > max depth seen in tests (33, Plummer N=5e4)
constexpr int DIRECT_TILE = 256;

template <typename Kernel>
int grid_for(Kernel kernel, int block, std::size_t n) {
    int device = 0, sms = 0, per_sm = 0;
    BH_CUDA_TRY(cudaGetDevice(&device));
    BH_CUDA_TRY(cudaDeviceGetAttribute(&sms, cudaDevAttrMultiProcessorCount, device));
    BH_CUDA_TRY(cudaOccupancyMaxActiveBlocksPerMultiprocessor(&per_sm, kernel, block, 0));
    const std::size_t needed = (n + block - 1) / block;
    const std::size_t resident = static_cast<std::size_t>(std::max(1, sms * per_sm));
    return static_cast<int>(std::max<std::size_t>(1, std::min(needed, resident)));
}

#define GRID_STRIDE(i, n) \
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < (n); i += blockDim.x * gridDim.x)

} // namespace

// ============================================================================
// Kernels
// ============================================================================

__global__ void __launch_bounds__(ELEM_BLOCK)
morton_kernel(int n, const float* __restrict__ px, const float* __restrict__ py,
              const float* __restrict__ pz, const Box* __restrict__ domain,
              std::uint64_t* __restrict__ keys, int* __restrict__ perm) {
    const Box box = *domain;
    GRID_STRIDE(i, n) {
        keys[i] = morton_key(px[i], py[i], pz[i], box);
        perm[i] = i;
    }
}

struct GatherArrays {
    const float* src[7];
    float* dst[7];
    const int* src_id;
    int* dst_id;
};

__global__ void __launch_bounds__(ELEM_BLOCK)
gather_kernel(int n, const int* __restrict__ perm, GatherArrays a) {
    GRID_STRIDE(i, n) {
        const int j = perm[i];
#pragma unroll
        for (int k = 0; k < 7; ++k) {
            a.dst[k][i] = a.src[k][j];
        }
        a.dst_id[i] = a.src_id[j];
    }
}

__global__ void __launch_bounds__(ELEM_BLOCK)
build_kernel(int n, const std::uint64_t* __restrict__ keys, int* __restrict__ left,
             int* __restrict__ right, int* __restrict__ parent) {
    GRID_STRIDE(i, n - 1) {
        build_internal_node(i, keys, n, left, right, parent);
    }
}

// TreeView arrays are deliberately not __restrict__/read-only: other threads
// write them and visibility relies on the acq_rel visit counters.
__global__ void __launch_bounds__(ELEM_BLOCK)
refit_kernel(int n, const float* __restrict__ px, const float* __restrict__ py,
             const float* __restrict__ pz, const float* __restrict__ mass, TreeView tree,
             float inv_theta2) {
    GRID_STRIDE(j, n) {
        refit_from_leaf(j, px, py, pz, mass, tree, inv_theta2);
    }
}

template <bool COLLECT_STATS>
__global__ void __launch_bounds__(FORCE_BLOCK)
force_kernel(int n, const int* __restrict__ left, const int* __restrict__ right,
             const float4* __restrict__ cm, const float4* __restrict__ lo,
             const float4* __restrict__ hi, float eps2, float G, float* __restrict__ fx,
             float* __restrict__ fy, float* __restrict__ fz,
             unsigned long long* __restrict__ counters) {
    TraversalCounts counts;
    GRID_STRIDE(i, n) {
        const float3 f = force_on<TRAVERSAL_STACK>(i, left, right, cm, lo, hi, n, eps2, G, counts);
        fx[i] = f.x;
        fy[i] = f.y;
        fz[i] = f.z;
    }
    if (counts.stack_overflows > 0) {
        atomicAdd(&counters[2], static_cast<unsigned long long>(counts.stack_overflows));
    }
    if constexpr (COLLECT_STATS) {
        using BlockReduce = cub::BlockReduce<unsigned long long, FORCE_BLOCK>;
        __shared__ typename BlockReduce::TempStorage storage;
        const unsigned long long pp = BlockReduce(storage).Sum(counts.particle_particle);
        __syncthreads();
        const unsigned long long pc = BlockReduce(storage).Sum(counts.particle_cell);
        if (threadIdx.x == 0) {
            atomicAdd(&counters[0], pp);
            atomicAdd(&counters[1], pc);
        }
    }
}

// Reference O(N^2): each block stages DIRECT_TILE sources in shared memory.
__global__ void __launch_bounds__(DIRECT_TILE)
direct_force_kernel(int n, const float* __restrict__ px, const float* __restrict__ py,
                    const float* __restrict__ pz, const float* __restrict__ mass, float eps2,
                    float G, float* __restrict__ fx, float* __restrict__ fy,
                    float* __restrict__ fz) {
    __shared__ float4 tile[DIRECT_TILE];
    const int tiles = (n + DIRECT_TILE - 1) / DIRECT_TILE;
    for (int base = blockIdx.x * blockDim.x; base < n; base += blockDim.x * gridDim.x) {
        const int i = base + threadIdx.x;
        const bool active = i < n;
        const float4 p = active ? float4{px[i], py[i], pz[i], mass[i]} : float4{0, 0, 0, 0};
        float3 f{0.0f, 0.0f, 0.0f};
        for (int t = 0; t < tiles; ++t) {
            const int j = t * DIRECT_TILE + threadIdx.x;
            tile[threadIdx.x] = j < n ? float4{px[j], py[j], pz[j], mass[j]} : float4{0, 0, 0, 0};
            __syncthreads();
            if (active) {
                const int count = min(DIRECT_TILE, n - t * DIRECT_TILE);
                for (int k = 0; k < count; ++k) {
                    if (t * DIRECT_TILE + k != i) {
                        accumulate(f, p, tile[k], eps2, G);
                    }
                }
            }
            __syncthreads();
        }
        if (active) {
            fx[i] = f.x;
            fy[i] = f.y;
            fz[i] = f.z;
        }
    }
}

__global__ void __launch_bounds__(ELEM_BLOCK)
kick_drift_kernel(int n, float* __restrict__ px, float* __restrict__ py, float* __restrict__ pz,
                  float* __restrict__ vx, float* __restrict__ vy, float* __restrict__ vz,
                  const float* __restrict__ fx, const float* __restrict__ fy,
                  const float* __restrict__ fz, const float* __restrict__ mass, float half_dt,
                  float dt) {
    GRID_STRIDE(i, n) {
        kick(i, vx, vy, vz, fx, fy, fz, mass, half_dt);
        drift(i, px, py, pz, vx, vy, vz, dt);
    }
}

__global__ void __launch_bounds__(ELEM_BLOCK)
kick_kernel(int n, float* __restrict__ vx, float* __restrict__ vy, float* __restrict__ vz,
            const float* __restrict__ fx, const float* __restrict__ fy,
            const float* __restrict__ fz, const float* __restrict__ mass, float half_dt) {
    GRID_STRIDE(i, n) {
        kick(i, vx, vy, vz, fx, fy, fz, mass, half_dt);
    }
}

// ============================================================================
// Host driver
// ============================================================================

bool BarnesHutGpu::device_available() noexcept {
    int count = 0;
    return cudaGetDeviceCount(&count) == cudaSuccess && count > 0;
}

BarnesHutGpu::BarnesHutGpu(const float* mass,
                           const float* pos_x, const float* pos_y, const float* pos_z,
                           const float* vel_x, const float* vel_y, const float* vel_z,
                           std::size_t count, const GpuConfig& config)
    : n_(count), config_(config) {
    if (count == 0 || count > static_cast<std::size_t>(INT32_MAX / 2)) {
        throw std::invalid_argument("BarnesHutGpu: particle count must be in [1, 2^30)");
    }
    const int n = static_cast<int>(n_);
    const std::size_t nodes = 2 * n_ - 1;
    const std::size_t internal = n_ > 1 ? n_ - 1 : 1;

    BH_CUDA_TRY(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
    for (auto& e : events_) {
        BH_CUDA_TRY(cudaEventCreate(&e));
    }

    px_ = device_alloc<float>(n_); py_ = device_alloc<float>(n_); pz_ = device_alloc<float>(n_);
    vx_ = device_alloc<float>(n_); vy_ = device_alloc<float>(n_); vz_ = device_alloc<float>(n_);
    m_ = device_alloc<float>(n_);
    fx_ = device_alloc<float>(n_); fy_ = device_alloc<float>(n_); fz_ = device_alloc<float>(n_);
    tmp_px_ = device_alloc<float>(n_); tmp_py_ = device_alloc<float>(n_); tmp_pz_ = device_alloc<float>(n_);
    tmp_vx_ = device_alloc<float>(n_); tmp_vy_ = device_alloc<float>(n_); tmp_vz_ = device_alloc<float>(n_);
    tmp_m_ = device_alloc<float>(n_);
    id_ = device_alloc<int>(n_); tmp_id_ = device_alloc<int>(n_);
    keys_in_ = device_alloc<std::uint64_t>(n_); keys_out_ = device_alloc<std::uint64_t>(n_);
    perm_in_ = device_alloc<int>(n_); perm_out_ = device_alloc<int>(n_);
    left_ = device_alloc<int>(internal); right_ = device_alloc<int>(internal);
    visits_ = device_alloc<int>(internal);
    parent_ = device_alloc<int>(nodes);
    cm_ = device_alloc<float4>(nodes); lo_ = device_alloc<float4>(nodes); hi_ = device_alloc<float4>(nodes);
    box_ = device_alloc<Box>(1);
    counters_ = device_alloc<unsigned long long>(3);

    // Temporary storage sized once for both CUB algorithms (graph-capturable)
    std::size_t reduce_bytes = 0, sort_bytes = 0;
    auto points = cuda::make_transform_iterator(cuda::counting_iterator<int>{0},
                                                PointBox{px_.get(), py_.get(), pz_.get()});
    BH_CUDA_TRY(cub::DeviceReduce::Reduce(nullptr, reduce_bytes, points, box_.get(), n, BoxMerge{},
                                          empty_box(), stream_));
    BH_CUDA_TRY(cub::DeviceRadixSort::SortPairs(nullptr, sort_bytes, keys_in_.get(), keys_out_.get(),
                                                perm_in_.get(), perm_out_.get(), n, 0, 63, stream_));
    temp_storage_bytes_ = std::max(reduce_bytes, sort_bytes);
    temp_storage_ = device_alloc<unsigned char>(temp_storage_bytes_);

    grid_elem_ = grid_for(refit_kernel, ELEM_BLOCK, n_);
    grid_force_ = grid_for(force_kernel<true>, FORCE_BLOCK, n_);

    // Upload; ids are the identity permutation of the caller's order
    const std::size_t bytes = n_ * sizeof(float);
    BH_CUDA_TRY(cudaMemcpyAsync(m_.get(), mass, bytes, cudaMemcpyHostToDevice, stream_));
    BH_CUDA_TRY(cudaMemcpyAsync(px_.get(), pos_x, bytes, cudaMemcpyHostToDevice, stream_));
    BH_CUDA_TRY(cudaMemcpyAsync(py_.get(), pos_y, bytes, cudaMemcpyHostToDevice, stream_));
    BH_CUDA_TRY(cudaMemcpyAsync(pz_.get(), pos_z, bytes, cudaMemcpyHostToDevice, stream_));
    BH_CUDA_TRY(cudaMemcpyAsync(vx_.get(), vel_x, bytes, cudaMemcpyHostToDevice, stream_));
    BH_CUDA_TRY(cudaMemcpyAsync(vy_.get(), vel_y, bytes, cudaMemcpyHostToDevice, stream_));
    BH_CUDA_TRY(cudaMemcpyAsync(vz_.get(), vel_z, bytes, cudaMemcpyHostToDevice, stream_));
    std::vector<int> ids(n_);
    for (int i = 0; i < n; ++i) ids[i] = i;
    BH_CUDA_TRY(cudaMemcpyAsync(id_.get(), ids.data(), n_ * sizeof(int), cudaMemcpyHostToDevice, stream_));

    // The root (node 0, or leaf 0 when n == 1) is the only node whose parent
    // is never written by build_kernel.
    BH_CUDA_TRY(cudaMemsetAsync(parent_.get(), 0xFF, nodes * sizeof(int), stream_));

    compute_forces();  // leapfrog needs a(x0)
    synchronize();     // ids vector must outlive the async copy
}

BarnesHutGpu::~BarnesHutGpu() {
    if (stream_) cudaStreamSynchronize(stream_);
    if (graph_exec_) cudaGraphExecDestroy(graph_exec_);
    if (graph_) cudaGraphDestroy(graph_);
    for (auto e : events_) {
        if (e) cudaEventDestroy(e);
    }
    if (stream_) cudaStreamDestroy(stream_);
}

void BarnesHutGpu::record(int slot) {
    BH_CUDA_TRY(cudaEventRecord(events_[slot], stream_));
}

void BarnesHutGpu::enqueue_tree() {
    const int n = static_cast<int>(n_);

    auto points = cuda::make_transform_iterator(cuda::counting_iterator<int>{0},
                                                PointBox{px_.get(), py_.get(), pz_.get()});
    std::size_t bytes = temp_storage_bytes_;
    BH_CUDA_TRY(cub::DeviceReduce::Reduce(temp_storage_.get(), bytes, points, box_.get(), n,
                                          BoxMerge{}, empty_box(), stream_));

    morton_kernel<<<grid_elem_, ELEM_BLOCK, 0, stream_>>>(n, px_.get(), py_.get(), pz_.get(), box_.get(),
                                                          keys_in_.get(), perm_in_.get());
    BH_CUDA_LAUNCH_CHECK();

    bytes = temp_storage_bytes_;
    BH_CUDA_TRY(cub::DeviceRadixSort::SortPairs(temp_storage_.get(), bytes, keys_in_.get(), keys_out_.get(),
                                                perm_in_.get(), perm_out_.get(), n, 0, 63, stream_));

    // Permute the particle state into key order (forces are recomputed)
    float* const state[7] = {px_.get(), py_.get(), pz_.get(), vx_.get(), vy_.get(), vz_.get(), m_.get()};
    float* const scratch[7] = {tmp_px_.get(), tmp_py_.get(), tmp_pz_.get(), tmp_vx_.get(), tmp_vy_.get(),
                               tmp_vz_.get(), tmp_m_.get()};
    GatherArrays g{};
    for (int k = 0; k < 7; ++k) {
        g.src[k] = state[k];
        g.dst[k] = scratch[k];
    }
    g.src_id = id_.get();
    g.dst_id = tmp_id_.get();
    gather_kernel<<<grid_elem_, ELEM_BLOCK, 0, stream_>>>(n, perm_out_.get(), g);
    BH_CUDA_LAUNCH_CHECK();
    // Copy back instead of swapping pointers so a captured graph stays valid
    const std::size_t fbytes = n_ * sizeof(float);
    for (int k = 0; k < 7; ++k) {
        BH_CUDA_TRY(cudaMemcpyAsync(state[k], scratch[k], fbytes, cudaMemcpyDeviceToDevice, stream_));
    }
    BH_CUDA_TRY(cudaMemcpyAsync(id_.get(), tmp_id_.get(), n_ * sizeof(int), cudaMemcpyDeviceToDevice, stream_));

    build_kernel<<<grid_elem_, ELEM_BLOCK, 0, stream_>>>(n, keys_out_.get(), left_.get(), right_.get(), parent_.get());
    BH_CUDA_LAUNCH_CHECK();

    if (n > 1) {
        BH_CUDA_TRY(cudaMemsetAsync(visits_.get(), 0, (n_ - 1) * sizeof(int), stream_));
    }
    const float inv_theta2 = config_.theta > 0.0f ? 1.0f / (config_.theta * config_.theta) : INFINITY;
    const TreeView tree{n, left_.get(), right_.get(), parent_.get(), cm_.get(), lo_.get(), hi_.get(), visits_.get()};
    refit_kernel<<<grid_elem_, ELEM_BLOCK, 0, stream_>>>(n, px_.get(), py_.get(), pz_.get(), m_.get(), tree, inv_theta2);
    BH_CUDA_LAUNCH_CHECK();
}

void BarnesHutGpu::enqueue_forces() {
    const int n = static_cast<int>(n_);
    BH_CUDA_TRY(cudaMemsetAsync(counters_.get(), 0, 3 * sizeof(unsigned long long), stream_));
    if (config_.collect_stats) {
        force_kernel<true><<<grid_force_, FORCE_BLOCK, 0, stream_>>>(
            n, left_.get(), right_.get(), cm_.get(), lo_.get(), hi_.get(), config_.epsilon_squared,
            config_.gravity, fx_.get(), fy_.get(), fz_.get(), counters_.get());
    } else {
        force_kernel<false><<<grid_force_, FORCE_BLOCK, 0, stream_>>>(
            n, left_.get(), right_.get(), cm_.get(), lo_.get(), hi_.get(), config_.epsilon_squared,
            config_.gravity, fx_.get(), fy_.get(), fz_.get(), counters_.get());
    }
    BH_CUDA_LAUNCH_CHECK();
}

void BarnesHutGpu::enqueue_step(bool timed) {
    const int n = static_cast<int>(n_);
    const float half_dt = 0.5f * config_.dt;

    if (timed) record(0);
    kick_drift_kernel<<<grid_elem_, ELEM_BLOCK, 0, stream_>>>(
        n, px_.get(), py_.get(), pz_.get(), vx_.get(), vy_.get(), vz_.get(), fx_.get(), fy_.get(),
        fz_.get(), m_.get(), half_dt, config_.dt);
    BH_CUDA_LAUNCH_CHECK();
    if (timed) record(1);
    enqueue_tree();
    if (timed) record(2);
    enqueue_forces();
    if (timed) record(3);
    kick_kernel<<<grid_elem_, ELEM_BLOCK, 0, stream_>>>(n, vx_.get(), vy_.get(), vz_.get(), fx_.get(),
                                                        fy_.get(), fz_.get(), m_.get(), half_dt);
    BH_CUDA_LAUNCH_CHECK();
    if (timed) record(4);
}

void BarnesHutGpu::capture_step_graph() {
    BH_CUDA_TRY(cudaStreamBeginCapture(stream_, cudaStreamCaptureModeThreadLocal));
    try {
        enqueue_step(false);
    } catch (...) {
        cudaGraph_t broken = nullptr;
        cudaStreamEndCapture(stream_, &broken);
        if (broken) cudaGraphDestroy(broken);
        throw;
    }
    BH_CUDA_TRY(cudaStreamEndCapture(stream_, &graph_));
    BH_CUDA_TRY(cudaGraphInstantiate(&graph_exec_, graph_, 0));
}

void BarnesHutGpu::compute_forces() {
    record(0);
    record(1);
    enqueue_tree();
    record(2);
    enqueue_forces();
    record(3);
    record(4);
    phase_timing_valid_ = true;
}

void BarnesHutGpu::compute_forces_direct() {
    const int n = static_cast<int>(n_);
    const int blocks = static_cast<int>((n_ + DIRECT_TILE - 1) / DIRECT_TILE);
    record(0);
    record(1);
    record(2);
    direct_force_kernel<<<blocks, DIRECT_TILE, 0, stream_>>>(n, px_.get(), py_.get(), pz_.get(), m_.get(),
                                                             config_.epsilon_squared, config_.gravity,
                                                             fx_.get(), fy_.get(), fz_.get());
    BH_CUDA_LAUNCH_CHECK();
    record(3);
    record(4);
    phase_timing_valid_ = true;
}

void BarnesHutGpu::step() {
    if (config_.use_graph) {
        if (!graph_exec_) {
            capture_step_graph();
        }
        record(0);
        BH_CUDA_TRY(cudaGraphLaunch(graph_exec_, stream_));
        record(4);
        phase_timing_valid_ = false;
    } else {
        enqueue_step(true);
        phase_timing_valid_ = true;
    }
}

GpuTiming BarnesHutGpu::synchronize() {
    BH_CUDA_TRY(cudaStreamSynchronize(stream_));
    GpuTiming t;
    BH_CUDA_TRY(cudaEventElapsedTime(&t.total_ms, events_[0], events_[4]));
    if (phase_timing_valid_) {
        float a = 0, b = 0;
        BH_CUDA_TRY(cudaEventElapsedTime(&a, events_[0], events_[1]));
        BH_CUDA_TRY(cudaEventElapsedTime(&t.tree_ms, events_[1], events_[2]));
        BH_CUDA_TRY(cudaEventElapsedTime(&t.force_ms, events_[2], events_[3]));
        BH_CUDA_TRY(cudaEventElapsedTime(&b, events_[3], events_[4]));
        t.integrate_ms = a + b;
    }
    return t;
}

void BarnesHutGpu::download(float* pos_x, float* pos_y, float* pos_z,
                            float* vel_x, float* vel_y, float* vel_z,
                            float* force_x, float* force_y, float* force_z) {
    std::vector<int> ids(n_);
    std::vector<float> buf(n_);
    BH_CUDA_TRY(cudaMemcpyAsync(ids.data(), id_.get(), n_ * sizeof(int), cudaMemcpyDeviceToHost, stream_));
    BH_CUDA_TRY(cudaStreamSynchronize(stream_));

    const std::pair<const float*, float*> fields[] = {
        {px_.get(), pos_x}, {py_.get(), pos_y}, {pz_.get(), pos_z},
        {vx_.get(), vel_x}, {vy_.get(), vel_y}, {vz_.get(), vel_z},
        {fx_.get(), force_x}, {fy_.get(), force_y}, {fz_.get(), force_z}};
    for (const auto& [src, dst] : fields) {
        if (!dst) continue;
        BH_CUDA_TRY(cudaMemcpyAsync(buf.data(), src, n_ * sizeof(float), cudaMemcpyDeviceToHost, stream_));
        BH_CUDA_TRY(cudaStreamSynchronize(stream_));
        for (std::size_t k = 0; k < n_; ++k) {
            dst[ids[k]] = buf[k];
        }
    }
}

GpuForceStats BarnesHutGpu::force_stats() {
    unsigned long long c[3] = {0, 0, 0};
    BH_CUDA_TRY(cudaMemcpyAsync(c, counters_.get(), sizeof(c), cudaMemcpyDeviceToHost, stream_));
    BH_CUDA_TRY(cudaStreamSynchronize(stream_));
    return GpuForceStats{c[0], c[1], c[2]};
}

} // namespace gpu
} // namespace barnes_hut
