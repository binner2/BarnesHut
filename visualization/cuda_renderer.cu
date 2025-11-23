#include "cuda_renderer.cuh"
#include <iostream>
#include <cmath>

namespace visualization {

// CUDA error checking macro
#define CUDA_CHECK(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(error) << std::endl; \
            return false; \
        } \
    } while(0)

#define CUDA_CHECK_VOID(call) \
    do { \
        cudaError_t error = call; \
        if (error != cudaSuccess) { \
            std::cerr << "CUDA error at " << __FILE__ << ":" << __LINE__ \
                      << " - " << cudaGetErrorString(error) << std::endl; \
        } \
    } while(0)

// ============================================================================
// CUDA Kernels
// ============================================================================

namespace kernels {

/**
 * @brief Convert particle data to vertex data with color mapping
 */
__global__ void prepare_vertices_kernel(
    const barnes_hut::Particle* particles,
    CudaRenderer::VertexData* vertices,
    size_t count,
    CudaRenderer::Config config
) {
    const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= count) return;

    const auto& particle = particles[idx];

    // Position
    vertices[idx].x = static_cast<float>(particle.position().x);
    vertices[idx].y = static_cast<float>(particle.position().y);
    vertices[idx].z = static_cast<float>(particle.position().z);

    // Color based on configuration
    float3 color;
    if (config.color_by_velocity) {
        // Calculate velocity magnitude
        const auto& vel = particle.velocity();
        float vel_mag = sqrtf(
            static_cast<float>(vel.x * vel.x + vel.y * vel.y + vel.z * vel.z)
        );
        color = color::velocity_to_color(vel_mag * config.velocity_scale, config.max_velocity);
    } else if (config.color_by_force) {
        // Calculate force magnitude
        const auto& force = particle.force();
        float force_mag = sqrtf(
            static_cast<float>(force.x * force.x + force.y * force.y + force.z * force.z)
        );
        color = color::force_to_color(force_mag * config.force_scale, config.max_force);
    } else {
        // Default white
        color = make_float3(1.0f, 1.0f, 1.0f);
    }

    vertices[idx].r = color.x;
    vertices[idx].g = color.y;
    vertices[idx].b = color.z;
    vertices[idx].a = 1.0f;
}

/**
 * @brief Compute velocity statistics using parallel reduction
 */
__global__ void compute_velocity_stats_kernel(
    const barnes_hut::Particle* particles,
    size_t count,
    float* min_vel,
    float* max_vel,
    float* sum_vel
) {
    __shared__ float shared_min[256];
    __shared__ float shared_max[256];
    __shared__ float shared_sum[256];

    const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t tid = threadIdx.x;

    // Initialize shared memory
    float local_min = INFINITY;
    float local_max = 0.0f;
    float local_sum = 0.0f;

    // Compute local values
    if (idx < count) {
        const auto& vel = particles[idx].velocity();
        float vel_mag = sqrtf(
            static_cast<float>(vel.x * vel.x + vel.y * vel.y + vel.z * vel.z)
        );
        local_min = vel_mag;
        local_max = vel_mag;
        local_sum = vel_mag;
    }

    shared_min[tid] = local_min;
    shared_max[tid] = local_max;
    shared_sum[tid] = local_sum;
    __syncthreads();

    // Parallel reduction in shared memory
    for (size_t s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            shared_min[tid] = fminf(shared_min[tid], shared_min[tid + s]);
            shared_max[tid] = fmaxf(shared_max[tid], shared_max[tid + s]);
            shared_sum[tid] += shared_sum[tid + s];
        }
        __syncthreads();
    }

    // Write block results to global memory
    if (tid == 0) {
        atomicMin((int*)min_vel, __float_as_int(shared_min[0]));
        atomicMax((int*)max_vel, __float_as_int(shared_max[0]));
        atomicAdd(sum_vel, shared_sum[0]);
    }
}

// Kernel launcher functions
void prepare_particle_vertices(
    const barnes_hut::Particle* particles,
    CudaRenderer::VertexData* vertices,
    size_t count,
    const CudaRenderer::Config& config,
    cudaStream_t stream
) {
    const int threads = 256;
    const int blocks = (count + threads - 1) / threads;

    prepare_vertices_kernel<<<blocks, threads, 0, stream>>>(
        particles, vertices, count, config
    );

    CUDA_CHECK_VOID(cudaGetLastError());
}

void compute_velocity_statistics(
    const barnes_hut::Particle* particles,
    size_t count,
    float* min_vel,
    float* max_vel,
    float* avg_vel,
    cudaStream_t stream
) {
    const int threads = 256;
    const int blocks = (count + threads - 1) / threads;

    // Initialize result buffers
    float init_min = INFINITY;
    float init_max = 0.0f;
    float init_sum = 0.0f;
    CUDA_CHECK_VOID(cudaMemcpyAsync(min_vel, &init_min, sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK_VOID(cudaMemcpyAsync(max_vel, &init_max, sizeof(float), cudaMemcpyHostToDevice, stream));
    CUDA_CHECK_VOID(cudaMemcpyAsync(avg_vel, &init_sum, sizeof(float), cudaMemcpyHostToDevice, stream));

    compute_velocity_stats_kernel<<<blocks, threads, 0, stream>>>(
        particles, count, min_vel, max_vel, avg_vel
    );

    CUDA_CHECK_VOID(cudaGetLastError());
}

} // namespace kernels

// ============================================================================
// CudaRenderer Implementation
// ============================================================================

CudaRenderer::CudaRenderer()
    : d_particles_(nullptr)
    , cuda_vbo_resource_(nullptr)
    , stream_(nullptr)
    , particle_count_(0)
    , initialized_(false)
{
}

CudaRenderer::~CudaRenderer() {
    cleanup();
}

bool CudaRenderer::initialize(size_t particle_count) {
    if (initialized_) {
        cleanup();
    }

    particle_count_ = particle_count;

    // Create CUDA stream for async operations
    CUDA_CHECK(cudaStreamCreate(&stream_));

    // Allocate device memory for particles
    if (!allocate_device_memory()) {
        cleanup();
        return false;
    }

    initialized_ = true;
    std::cout << "CUDA Renderer initialized for " << particle_count << " particles\n";

    return true;
}

bool CudaRenderer::register_gl_buffer(unsigned int vbo) {
    if (!initialized_) {
        std::cerr << "CudaRenderer not initialized\n";
        return false;
    }

    // Unregister previous buffer if exists
    if (cuda_vbo_resource_) {
        unregister_gl_buffer();
    }

    // Register OpenGL buffer with CUDA
    CUDA_CHECK(cudaGraphicsGLRegisterBuffer(
        &cuda_vbo_resource_,
        vbo,
        cudaGraphicsMapFlagsWriteDiscard
    ));

    std::cout << "OpenGL VBO registered with CUDA\n";
    return true;
}

void CudaRenderer::unregister_gl_buffer() {
    if (cuda_vbo_resource_) {
        CUDA_CHECK_VOID(cudaGraphicsUnregisterResource(cuda_vbo_resource_));
        cuda_vbo_resource_ = nullptr;
    }
}

bool CudaRenderer::update_particles(
    const barnes_hut::Particle* particles,
    size_t count,
    const Config& config
) {
    if (!initialized_ || !cuda_vbo_resource_) {
        return false;
    }

    cudaEvent_t start, stop;
    CUDA_CHECK(cudaEventCreate(&start));
    CUDA_CHECK(cudaEventCreate(&stop));
    CUDA_CHECK(cudaEventRecord(start, stream_));

    config_ = config;

    // Copy particle data to device
    CUDA_CHECK(cudaMemcpyAsync(
        d_particles_,
        particles,
        count * sizeof(barnes_hut::Particle),
        cudaMemcpyHostToDevice,
        stream_
    ));

    // Map OpenGL buffer for CUDA access
    CUDA_CHECK(cudaGraphicsMapResources(1, &cuda_vbo_resource_, stream_));

    VertexData* d_vertices = nullptr;
    size_t num_bytes = 0;
    CUDA_CHECK(cudaGraphicsResourceGetMappedPointer(
        (void**)&d_vertices,
        &num_bytes,
        cuda_vbo_resource_
    ));

    // Launch kernel to prepare vertex data
    kernels::prepare_particle_vertices(
        static_cast<const barnes_hut::Particle*>(d_particles_),
        d_vertices,
        count,
        config_,
        stream_
    );

    // Unmap OpenGL buffer
    CUDA_CHECK(cudaGraphicsUnmapResources(1, &cuda_vbo_resource_, stream_));

    // Compute statistics (optional, can be disabled for performance)
    if constexpr (false) {  // Disable for now
        float* d_min, * d_max, * d_avg;
        CUDA_CHECK(cudaMallocAsync(&d_min, sizeof(float), stream_));
        CUDA_CHECK(cudaMallocAsync(&d_max, sizeof(float), stream_));
        CUDA_CHECK(cudaMallocAsync(&d_avg, sizeof(float), stream_));

        kernels::compute_velocity_statistics(
            static_cast<const barnes_hut::Particle*>(d_particles_),
            count, d_min, d_max, d_avg, stream_
        );

        CUDA_CHECK(cudaMemcpyAsync(&stats_.min_velocity, d_min, sizeof(float), cudaMemcpyDeviceToHost, stream_));
        CUDA_CHECK(cudaMemcpyAsync(&stats_.max_velocity, d_max, sizeof(float), cudaMemcpyDeviceToHost, stream_));
        CUDA_CHECK(cudaMemcpyAsync(&stats_.avg_velocity, d_avg, sizeof(float), cudaMemcpyDeviceToHost, stream_));

        CUDA_CHECK(cudaFreeAsync(d_min, stream_));
        CUDA_CHECK(cudaFreeAsync(d_max, stream_));
        CUDA_CHECK(cudaFreeAsync(d_avg, stream_));
    }

    CUDA_CHECK(cudaEventRecord(stop, stream_));
    CUDA_CHECK(cudaEventSynchronize(stop));

    float milliseconds = 0;
    CUDA_CHECK(cudaEventElapsedTime(&milliseconds, start, stop));
    stats_.last_update_time_ms = milliseconds;
    stats_.particle_count = count;

    CUDA_CHECK(cudaEventDestroy(start));
    CUDA_CHECK(cudaEventDestroy(stop));

    return true;
}

void CudaRenderer::cleanup() {
    if (cuda_vbo_resource_) {
        unregister_gl_buffer();
    }

    free_device_memory();

    if (stream_) {
        CUDA_CHECK_VOID(cudaStreamDestroy(stream_));
        stream_ = nullptr;
    }

    initialized_ = false;
}

bool CudaRenderer::allocate_device_memory() {
    CUDA_CHECK(cudaMalloc(
        &d_particles_,
        particle_count_ * sizeof(barnes_hut::Particle)
    ));

    return true;
}

void CudaRenderer::free_device_memory() {
    if (d_particles_) {
        CUDA_CHECK_VOID(cudaFree(d_particles_));
        d_particles_ = nullptr;
    }
}

} // namespace visualization
