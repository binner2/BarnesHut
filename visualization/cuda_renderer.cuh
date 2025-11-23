#pragma once

#include <cuda_runtime.h>
#include <cuda_gl_interop.h>
#include "../particle.h"
#include "../vektor.h"

namespace visualization {

/**
 * @brief CUDA-based particle renderer with OpenGL interoperability
 *
 * This class provides GPU-accelerated particle rendering using CUDA kernels
 * with zero-copy data transfer via CUDA-OpenGL interop.
 */
class CudaRenderer {
public:
    /**
     * @brief Vertex data structure for rendering (position + color)
     */
    struct alignas(16) VertexData {
        float x, y, z;      // Position
        float r, g, b, a;   // Color with alpha
    };

    /**
     * @brief Rendering configuration
     */
    struct Config {
        bool color_by_velocity = true;   // Color particles by velocity
        bool color_by_force = false;      // Color particles by force magnitude
        float point_size = 2.0f;          // Base point size
        float velocity_scale = 1.0f;      // Velocity color scale
        float force_scale = 1.0f;         // Force color scale
        float max_velocity = 100.0f;      // Maximum velocity for color mapping
        float max_force = 10000.0f;       // Maximum force for color mapping
    };

    CudaRenderer();
    ~CudaRenderer();

    // Delete copy operations
    CudaRenderer(const CudaRenderer&) = delete;
    CudaRenderer& operator=(const CudaRenderer&) = delete;

    /**
     * @brief Initialize CUDA renderer with particle count
     * @param particle_count Number of particles to render
     * @return true if initialization successful
     */
    bool initialize(size_t particle_count);

    /**
     * @brief Register OpenGL buffer for CUDA interop
     * @param vbo OpenGL Vertex Buffer Object ID
     * @return true if registration successful
     */
    bool register_gl_buffer(unsigned int vbo);

    /**
     * @brief Unregister OpenGL buffer
     */
    void unregister_gl_buffer();

    /**
     * @brief Update particle data on GPU and prepare for rendering
     * @param particles Span of particles to render
     * @param config Rendering configuration
     * @return true if update successful
     */
    bool update_particles(
        const barnes_hut::Particle* particles,
        size_t count,
        const Config& config
    );

    /**
     * @brief Get rendering configuration (mutable)
     */
    Config& get_config() { return config_; }

    /**
     * @brief Get rendering statistics
     */
    struct Statistics {
        double last_update_time_ms = 0.0;
        size_t particle_count = 0;
        float min_velocity = 0.0f;
        float max_velocity = 0.0f;
        float avg_velocity = 0.0f;
    };

    const Statistics& get_statistics() const { return stats_; }

    /**
     * @brief Cleanup resources
     */
    void cleanup();

private:
    // CUDA resources
    void* d_particles_;                    // Device particle data
    cudaGraphicsResource* cuda_vbo_resource_; // CUDA-GL interop resource
    cudaStream_t stream_;                  // CUDA stream for async operations

    // Configuration
    Config config_;
    Statistics stats_;
    size_t particle_count_;
    bool initialized_;

    // Helper methods
    bool allocate_device_memory();
    void free_device_memory();
};

// CUDA kernel declarations (implemented in .cu file)
namespace kernels {

/**
 * @brief Convert particle data to vertex data with color mapping
 * @param particles Input particle array
 * @param vertices Output vertex array (OpenGL buffer)
 * @param count Number of particles
 * @param config Rendering configuration
 */
void prepare_particle_vertices(
    const barnes_hut::Particle* particles,
    CudaRenderer::VertexData* vertices,
    size_t count,
    const CudaRenderer::Config& config,
    cudaStream_t stream
);

/**
 * @brief Compute velocity statistics on GPU
 */
void compute_velocity_statistics(
    const barnes_hut::Particle* particles,
    size_t count,
    float* min_vel,
    float* max_vel,
    float* avg_vel,
    cudaStream_t stream
);

} // namespace kernels

// Utility functions for color mapping
namespace color {

/**
 * @brief Map velocity to color (blue -> green -> yellow -> red)
 * @param velocity Velocity magnitude [0, max_velocity]
 * @param max_velocity Maximum velocity for normalization
 * @return RGB color (0-1 range)
 */
__device__ __host__ inline float3 velocity_to_color(float velocity, float max_velocity) {
    float t = fminf(velocity / max_velocity, 1.0f);

    float3 color;
    if (t < 0.25f) {
        // Blue to Cyan
        float s = t * 4.0f;
        color = make_float3(0.0f, s, 1.0f);
    } else if (t < 0.5f) {
        // Cyan to Green
        float s = (t - 0.25f) * 4.0f;
        color = make_float3(0.0f, 1.0f, 1.0f - s);
    } else if (t < 0.75f) {
        // Green to Yellow
        float s = (t - 0.5f) * 4.0f;
        color = make_float3(s, 1.0f, 0.0f);
    } else {
        // Yellow to Red
        float s = (t - 0.75f) * 4.0f;
        color = make_float3(1.0f, 1.0f - s, 0.0f);
    }

    return color;
}

/**
 * @brief Map force magnitude to color (white -> yellow -> orange -> red -> magenta)
 */
__device__ __host__ inline float3 force_to_color(float force, float max_force) {
    float t = fminf(force / max_force, 1.0f);

    float3 color;
    if (t < 0.33f) {
        // White to Yellow
        float s = t * 3.0f;
        color = make_float3(1.0f, 1.0f, 1.0f - s);
    } else if (t < 0.66f) {
        // Yellow to Red
        float s = (t - 0.33f) * 3.0f;
        color = make_float3(1.0f, 1.0f - s, 0.0f);
    } else {
        // Red to Magenta
        float s = (t - 0.66f) * 3.0f;
        color = make_float3(1.0f, 0.0f, s);
    }

    return color;
}

} // namespace color

} // namespace visualization
