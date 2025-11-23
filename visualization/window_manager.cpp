#include "window_manager.h"
#include <glad/glad.h>
#include <iostream>

namespace visualization {

WindowManager::WindowManager()
    : window_(nullptr)
    , last_mouse_x_(0.0)
    , last_mouse_y_(0.0)
    , left_button_pressed_(false)
    , right_button_pressed_(false)
    , first_mouse_(true)
    , last_frame_time_(0.0)
    , delta_time_(0.0)
    , fps_(0.0)
    , frame_count_(0)
    , fps_update_time_(0.0)
    , fps_frame_count_(0)
    , initialized_(false)
{
}

WindowManager::~WindowManager() {
    cleanup();
}

bool WindowManager::initialize(const Config& config) {
    if (initialized_) {
        cleanup();
    }

    config_ = config;

    // Set GLFW error callback
    glfwSetErrorCallback(error_callback);

    // Initialize GLFW
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return false;
    }

    // Configure GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, config.msaa_samples);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // Create window
    GLFWmonitor* monitor = config.fullscreen ? glfwGetPrimaryMonitor() : nullptr;
    window_ = glfwCreateWindow(config.width, config.height, config.title.c_str(), monitor, nullptr);

    if (!window_) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return false;
    }

    // Make context current
    glfwMakeContextCurrent(window_);

    // Enable/disable VSync
    glfwSwapInterval(config.vsync ? 1 : 0);

    // Load OpenGL functions with GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD\n";
        glfwDestroyWindow(window_);
        glfwTerminate();
        return false;
    }

    // Set user pointer for callbacks
    glfwSetWindowUserPointer(window_, this);

    // Register callbacks
    glfwSetCursorPosCallback(window_, cursor_position_callback);
    glfwSetMouseButtonCallback(window_, mouse_button_callback);
    glfwSetScrollCallback(window_, scroll_callback);
    glfwSetKeyCallback(window_, key_callback);
    glfwSetFramebufferSizeCallback(window_, framebuffer_size_callback);

    // Initialize timing
    last_frame_time_ = glfwGetTime();
    fps_update_time_ = last_frame_time_;

    initialized_ = true;

    std::cout << "Window initialized: " << config.width << "x" << config.height << "\n";

    return true;
}

bool WindowManager::should_close() const {
    return window_ && glfwWindowShouldClose(window_);
}

void WindowManager::begin_frame() {
    // Update timing
    double current_time = glfwGetTime();
    delta_time_ = current_time - last_frame_time_;
    last_frame_time_ = current_time;

    // Update FPS counter (every second)
    fps_frame_count_++;
    if (current_time - fps_update_time_ >= 1.0) {
        fps_ = fps_frame_count_ / (current_time - fps_update_time_);
        fps_update_time_ = current_time;
        fps_frame_count_ = 0;

        // Update window title with FPS
        std::string title = config_.title + " | FPS: " + std::to_string(static_cast<int>(fps_));
        glfwSetWindowTitle(window_, title.c_str());
    }

    frame_count_++;
}

void WindowManager::end_frame() {
    // Swap buffers and poll events
    glfwSwapBuffers(window_);
    glfwPollEvents();
}

void WindowManager::get_window_size(int& width, int& height) const {
    glfwGetWindowSize(window_, &width, &height);
}

void WindowManager::get_framebuffer_size(int& width, int& height) const {
    glfwGetFramebufferSize(window_, &width, &height);
}

void WindowManager::set_input_callbacks(const InputCallbacks& callbacks) {
    callbacks_ = callbacks;
}

void WindowManager::cleanup() {
    if (window_) {
        glfwDestroyWindow(window_);
        window_ = nullptr;
    }

    if (initialized_) {
        glfwTerminate();
        initialized_ = false;
    }
}

// ============================================================================
// Static GLFW Callbacks
// ============================================================================

void WindowManager::cursor_position_callback(GLFWwindow* window, double xpos, double ypos) {
    auto* manager = static_cast<WindowManager*>(glfwGetWindowUserPointer(window));
    if (manager) {
        manager->handle_cursor_position(xpos, ypos);
    }
}

void WindowManager::mouse_button_callback(GLFWwindow* window, int button, int action, int mods) {
    auto* manager = static_cast<WindowManager*>(glfwGetWindowUserPointer(window));
    if (manager) {
        manager->handle_mouse_button(button, action, mods);
    }
}

void WindowManager::scroll_callback(GLFWwindow* window, double xoffset, double yoffset) {
    auto* manager = static_cast<WindowManager*>(glfwGetWindowUserPointer(window));
    if (manager) {
        manager->handle_scroll(xoffset, yoffset);
    }
}

void WindowManager::key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    auto* manager = static_cast<WindowManager*>(glfwGetWindowUserPointer(window));
    if (manager) {
        manager->handle_key(key, scancode, action, mods);
    }
}

void WindowManager::framebuffer_size_callback(GLFWwindow* window, int width, int height) {
    auto* manager = static_cast<WindowManager*>(glfwGetWindowUserPointer(window));
    if (manager) {
        manager->handle_framebuffer_resize(width, height);
    }
}

void WindowManager::error_callback(int error, const char* description) {
    std::cerr << "GLFW Error " << error << ": " << description << "\n";
}

// ============================================================================
// Instance Callback Handlers
// ============================================================================

void WindowManager::handle_cursor_position(double xpos, double ypos) {
    if (first_mouse_) {
        last_mouse_x_ = xpos;
        last_mouse_y_ = ypos;
        first_mouse_ = false;
        return;
    }

    float delta_x = static_cast<float>(xpos - last_mouse_x_);
    float delta_y = static_cast<float>(ypos - last_mouse_y_);

    last_mouse_x_ = xpos;
    last_mouse_y_ = ypos;

    if (callbacks_.on_mouse_move) {
        callbacks_.on_mouse_move(delta_x, delta_y, left_button_pressed_, right_button_pressed_);
    }
}

void WindowManager::handle_mouse_button(int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        left_button_pressed_ = (action == GLFW_PRESS);
    } else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
        right_button_pressed_ = (action == GLFW_PRESS);
    }
}

void WindowManager::handle_scroll(double xoffset, double yoffset) {
    if (callbacks_.on_mouse_wheel) {
        callbacks_.on_mouse_wheel(static_cast<float>(yoffset));
    }
}

void WindowManager::handle_key(int key, int scancode, int action, int mods) {
    // Handle common keys
    if (action == GLFW_PRESS) {
        if (key == GLFW_KEY_ESCAPE) {
            glfwSetWindowShouldClose(window_, GLFW_TRUE);
        } else if (key == GLFW_KEY_F11) {
            // Toggle fullscreen
            config_.fullscreen = !config_.fullscreen;
            // Implementation for fullscreen toggle can be added here
        }

        if (callbacks_.on_key_press) {
            callbacks_.on_key_press(key);
        }
    }
}

void WindowManager::handle_framebuffer_resize(int width, int height) {
    glViewport(0, 0, width, height);
}

} // namespace visualization
