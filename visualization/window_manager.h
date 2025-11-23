#pragma once

#include <GLFW/glfw3.h>
#include <string>
#include <functional>

namespace visualization {

/**
 * @brief Window manager for GLFW-based visualization
 *
 * Handles window creation, input events, and frame timing.
 */
class WindowManager {
public:
    /**
     * @brief Window configuration
     */
    struct Config {
        int width = 1920;
        int height = 1080;
        std::string title = "Barnes-Hut N-Body Simulation - CUDA Visualization";
        bool fullscreen = false;
        bool vsync = true;
        int msaa_samples = 4;  // Multisample anti-aliasing
    };

    /**
     * @brief Input callbacks
     */
    struct InputCallbacks {
        std::function<void(float, float, bool, bool)> on_mouse_move;  // delta_x, delta_y, left_btn, right_btn
        std::function<void(float)> on_mouse_wheel;                    // delta
        std::function<void(int)> on_key_press;                        // key code
    };

    WindowManager();
    ~WindowManager();

    // Delete copy operations
    WindowManager(const WindowManager&) = delete;
    WindowManager& operator=(const WindowManager&) = delete;

    /**
     * @brief Initialize window and OpenGL context
     * @param config Window configuration
     * @return true if initialization successful
     */
    bool initialize(const Config& config);

    /**
     * @brief Check if window should close
     */
    bool should_close() const;

    /**
     * @brief Begin frame (clear and prepare for rendering)
     */
    void begin_frame();

    /**
     * @brief End frame (swap buffers and poll events)
     */
    void end_frame();

    /**
     * @brief Get window dimensions
     */
    void get_window_size(int& width, int& height) const;

    /**
     * @brief Get framebuffer dimensions (may differ from window size on high-DPI displays)
     */
    void get_framebuffer_size(int& width, int& height) const;

    /**
     * @brief Set input callbacks
     */
    void set_input_callbacks(const InputCallbacks& callbacks);

    /**
     * @brief Get frame timing information
     */
    double get_delta_time() const { return delta_time_; }
    double get_fps() const { return fps_; }
    uint64_t get_frame_count() const { return frame_count_; }

    /**
     * @brief Get GLFW window handle
     */
    GLFWwindow* get_handle() { return window_; }

    /**
     * @brief Cleanup resources
     */
    void cleanup();

private:
    // GLFW callbacks (static, redirect to instance methods)
    static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos);
    static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods);
    static void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
    static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods);
    static void framebuffer_size_callback(GLFWwindow* window, int width, int height);
    static void error_callback(int error, const char* description);

    // Instance callback handlers
    void handle_cursor_position(double xpos, double ypos);
    void handle_mouse_button(int button, int action, int mods);
    void handle_scroll(double xoffset, double yoffset);
    void handle_key(int key, int scancode, int action, int mods);
    void handle_framebuffer_resize(int width, int height);

    GLFWwindow* window_;
    Config config_;
    InputCallbacks callbacks_;

    // Mouse state
    double last_mouse_x_;
    double last_mouse_y_;
    bool left_button_pressed_;
    bool right_button_pressed_;
    bool first_mouse_;

    // Timing
    double last_frame_time_;
    double delta_time_;
    double fps_;
    uint64_t frame_count_;
    double fps_update_time_;
    uint64_t fps_frame_count_;

    bool initialized_;
};

} // namespace visualization
