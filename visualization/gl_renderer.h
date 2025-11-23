#pragma once

#include <string>
#include <memory>
#include <glad/glad.h>
#include <glm/glm.hpp>

namespace visualization {

/**
 * @brief Modern OpenGL renderer for particle systems
 *
 * Uses OpenGL 4.5+ features for efficient particle rendering with
 * CUDA interoperability support.
 */
class GLRenderer {
public:
    /**
     * @brief Rendering configuration
     */
    struct Config {
        float point_size = 3.0f;
        float point_scale = 100.0f;
        float brightness = 1.2f;
        float alpha = 0.9f;
        bool enable_glow = true;
        bool enable_blend = true;
        bool enable_depth_test = true;
        glm::vec3 clear_color = glm::vec3(0.02f, 0.02f, 0.05f);  // Dark blue background
    };

    GLRenderer();
    ~GLRenderer();

    // Delete copy operations
    GLRenderer(const GLRenderer&) = delete;
    GLRenderer& operator=(const GLRenderer&) = delete;

    /**
     * @brief Initialize OpenGL renderer
     * @param particle_count Number of particles to render
     * @return true if initialization successful
     */
    bool initialize(size_t particle_count);

    /**
     * @brief Set camera matrices
     * @param view View matrix
     * @param projection Projection matrix
     * @param camera_position Camera position in world space
     */
    void set_camera(const glm::mat4& view, const glm::mat4& projection, const glm::vec3& camera_position);

    /**
     * @brief Render particles
     * @param particle_count Number of particles to render
     */
    void render(size_t particle_count);

    /**
     * @brief Get OpenGL VBO handle (for CUDA interop)
     */
    GLuint get_vbo() const { return vbo_; }

    /**
     * @brief Get/set rendering configuration
     */
    Config& get_config() { return config_; }
    const Config& get_config() const { return config_; }

    /**
     * @brief Cleanup resources
     */
    void cleanup();

private:
    // Shader management
    bool load_shaders();
    bool compile_shader(GLuint shader, const std::string& source);
    bool link_program(GLuint program);
    std::string read_shader_file(const std::string& filepath);

    // OpenGL resources
    GLuint vao_;           // Vertex Array Object
    GLuint vbo_;           // Vertex Buffer Object
    GLuint shader_program_; // Shader program

    // Shader uniform locations
    GLint u_view_projection_;
    GLint u_camera_position_;
    GLint u_point_size_;
    GLint u_point_scale_;
    GLint u_brightness_;
    GLint u_enable_glow_;
    GLint u_alpha_;

    // Configuration
    Config config_;
    size_t particle_count_;
    bool initialized_;

    // Camera state
    glm::mat4 view_projection_;
    glm::vec3 camera_position_;
};

} // namespace visualization
