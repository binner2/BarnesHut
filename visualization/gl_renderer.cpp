#include "gl_renderer.h"
#include "cuda_renderer.cuh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <glm/gtc/type_ptr.hpp>

namespace visualization {

GLRenderer::GLRenderer()
    : vao_(0)
    , vbo_(0)
    , shader_program_(0)
    , particle_count_(0)
    , initialized_(false)
    , camera_position_(0.0f)
{
}

GLRenderer::~GLRenderer() {
    cleanup();
}

bool GLRenderer::initialize(size_t particle_count) {
    if (initialized_) {
        cleanup();
    }

    particle_count_ = particle_count;

    // Enable OpenGL debug output (if available)
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);

    // Create and bind VAO
    glGenVertexArrays(1, &vao_);
    glBindVertexArray(vao_);

    // Create VBO
    glGenBuffers(1, &vbo_);
    glBindBuffer(GL_ARRAY_BUFFER, vbo_);

    // Allocate buffer storage (immutable storage for better performance)
    const size_t buffer_size = particle_count * sizeof(CudaRenderer::VertexData);
    glBufferData(GL_ARRAY_BUFFER, buffer_size, nullptr, GL_DYNAMIC_DRAW);

    // Setup vertex attributes
    // Position (location = 0): vec3
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE,
                         sizeof(CudaRenderer::VertexData),
                         (void*)offsetof(CudaRenderer::VertexData, x));
    glEnableVertexAttribArray(0);

    // Color (location = 1): vec4
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE,
                         sizeof(CudaRenderer::VertexData),
                         (void*)offsetof(CudaRenderer::VertexData, r));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    // Load and compile shaders
    if (!load_shaders()) {
        std::cerr << "Failed to load shaders\n";
        cleanup();
        return false;
    }

    // Get uniform locations
    u_view_projection_ = glGetUniformLocation(shader_program_, "u_view_projection");
    u_camera_position_ = glGetUniformLocation(shader_program_, "u_camera_position");
    u_point_size_ = glGetUniformLocation(shader_program_, "u_point_size");
    u_point_scale_ = glGetUniformLocation(shader_program_, "u_point_scale");
    u_brightness_ = glGetUniformLocation(shader_program_, "u_brightness");
    u_enable_glow_ = glGetUniformLocation(shader_program_, "u_enable_glow");
    u_alpha_ = glGetUniformLocation(shader_program_, "u_alpha");

    // Configure OpenGL state
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_POINT_SPRITE);

    if (config_.enable_blend) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        // For additive blending (glow effect):
        // glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    }

    if (config_.enable_depth_test) {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
    }

    initialized_ = true;
    std::cout << "OpenGL Renderer initialized for " << particle_count << " particles\n";
    std::cout << "  OpenGL Version: " << glGetString(GL_VERSION) << "\n";
    std::cout << "  GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << "\n";
    std::cout << "  Renderer: " << glGetString(GL_RENDERER) << "\n";

    return true;
}

void GLRenderer::set_camera(const glm::mat4& view, const glm::mat4& projection, const glm::vec3& camera_position) {
    view_projection_ = projection * view;
    camera_position_ = camera_position;
}

void GLRenderer::render(size_t particle_count) {
    if (!initialized_) {
        return;
    }

    // Clear buffers
    glClearColor(config_.clear_color.r, config_.clear_color.g, config_.clear_color.b, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Use shader program
    glUseProgram(shader_program_);

    // Update uniforms
    glUniformMatrix4fv(u_view_projection_, 1, GL_FALSE, glm::value_ptr(view_projection_));
    glUniform3fv(u_camera_position_, 1, glm::value_ptr(camera_position_));
    glUniform1f(u_point_size_, config_.point_size);
    glUniform1f(u_point_scale_, config_.point_scale);
    glUniform1f(u_brightness_, config_.brightness);
    glUniform1i(u_enable_glow_, config_.enable_glow ? 1 : 0);
    glUniform1f(u_alpha_, config_.alpha);

    // Bind VAO and draw
    glBindVertexArray(vao_);
    glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(particle_count));
    glBindVertexArray(0);

    // Check for errors
    GLenum error = glGetError();
    if (error != GL_NO_ERROR) {
        std::cerr << "OpenGL error: " << error << "\n";
    }
}

void GLRenderer::cleanup() {
    if (shader_program_) {
        glDeleteProgram(shader_program_);
        shader_program_ = 0;
    }

    if (vbo_) {
        glDeleteBuffers(1, &vbo_);
        vbo_ = 0;
    }

    if (vao_) {
        glDeleteVertexArrays(1, &vao_);
        vao_ = 0;
    }

    initialized_ = false;
}

bool GLRenderer::load_shaders() {
    // Read shader source files
    std::string vertex_source = read_shader_file("shaders/particle.vert");
    std::string fragment_source = read_shader_file("shaders/particle.frag");

    if (vertex_source.empty() || fragment_source.empty()) {
        std::cerr << "Failed to read shader files\n";
        return false;
    }

    // Create and compile vertex shader
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    if (!compile_shader(vertex_shader, vertex_source)) {
        glDeleteShader(vertex_shader);
        return false;
    }

    // Create and compile fragment shader
    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    if (!compile_shader(fragment_shader, fragment_source)) {
        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
        return false;
    }

    // Create and link shader program
    shader_program_ = glCreateProgram();
    glAttachShader(shader_program_, vertex_shader);
    glAttachShader(shader_program_, fragment_shader);

    if (!link_program(shader_program_)) {
        glDeleteShader(vertex_shader);
        glDeleteShader(fragment_shader);
        glDeleteProgram(shader_program_);
        shader_program_ = 0;
        return false;
    }

    // Shaders can be deleted after linking
    glDeleteShader(vertex_shader);
    glDeleteShader(fragment_shader);

    std::cout << "Shaders compiled and linked successfully\n";
    return true;
}

bool GLRenderer::compile_shader(GLuint shader, const std::string& source) {
    const char* source_cstr = source.c_str();
    glShaderSource(shader, 1, &source_cstr, nullptr);
    glCompileShader(shader);

    GLint success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);

    if (!success) {
        char info_log[1024];
        glGetShaderInfoLog(shader, sizeof(info_log), nullptr, info_log);
        std::cerr << "Shader compilation error:\n" << info_log << "\n";
        return false;
    }

    return true;
}

bool GLRenderer::link_program(GLuint program) {
    glLinkProgram(program);

    GLint success;
    glGetProgramiv(program, GL_LINK_STATUS, &success);

    if (!success) {
        char info_log[1024];
        glGetProgramInfoLog(program, sizeof(info_log), nullptr, info_log);
        std::cerr << "Shader linking error:\n" << info_log << "\n";
        return false;
    }

    return true;
}

std::string GLRenderer::read_shader_file(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open shader file: " << filepath << "\n";
        return "";
    }

    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

} // namespace visualization
