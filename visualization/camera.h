#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace visualization {

/**
 * @brief Orbit camera for 3D visualization
 *
 * Provides intuitive camera controls:
 * - Mouse drag: Orbit around target
 * - Mouse wheel: Zoom in/out
 * - Right mouse + drag: Pan
 */
class Camera {
public:
    Camera();

    /**
     * @brief Initialize camera with target position
     * @param target Center point to orbit around
     * @param distance Initial distance from target
     */
    void initialize(const glm::vec3& target, float distance = 50.0f);

    /**
     * @brief Handle mouse movement for camera control
     * @param delta_x Mouse movement in X direction
     * @param delta_y Mouse movement in Y direction
     * @param left_button Left mouse button pressed
     * @param right_button Right mouse button pressed
     */
    void handle_mouse_move(float delta_x, float delta_y, bool left_button, bool right_button);

    /**
     * @brief Handle mouse wheel for zoom
     * @param delta Wheel delta (positive = zoom in)
     */
    void handle_mouse_wheel(float delta);

    /**
     * @brief Update camera matrices
     */
    void update();

    /**
     * @brief Set viewport size (for projection matrix)
     * @param width Viewport width in pixels
     * @param height Viewport height in pixels
     */
    void set_viewport(int width, int height);

    /**
     * @brief Get view matrix
     */
    const glm::mat4& get_view_matrix() const { return view_matrix_; }

    /**
     * @brief Get projection matrix
     */
    const glm::mat4& get_projection_matrix() const { return projection_matrix_; }

    /**
     * @brief Get camera position in world space
     */
    glm::vec3 get_position() const { return position_; }

    /**
     * @brief Get target position
     */
    const glm::vec3& get_target() const { return target_; }

    /**
     * @brief Reset camera to initial position
     */
    void reset();

    /**
     * @brief Camera configuration
     */
    struct Config {
        float fov = 45.0f;              // Field of view in degrees
        float near_plane = 0.1f;        // Near clipping plane
        float far_plane = 1000.0f;      // Far clipping plane
        float rotation_speed = 0.005f;  // Mouse rotation sensitivity
        float pan_speed = 0.05f;        // Mouse pan sensitivity
        float zoom_speed = 2.0f;        // Mouse wheel sensitivity
        float min_distance = 1.0f;      // Minimum zoom distance
        float max_distance = 500.0f;    // Maximum zoom distance
    };

    Config& get_config() { return config_; }
    const Config& get_config() const { return config_; }

private:
    void update_position();
    void update_matrices();

    // Camera state
    glm::vec3 target_;          // Point camera orbits around
    glm::vec3 position_;        // Current camera position
    float distance_;            // Distance from target
    float azimuth_;             // Horizontal angle (radians)
    float elevation_;           // Vertical angle (radians)

    // Matrices
    glm::mat4 view_matrix_;
    glm::mat4 projection_matrix_;

    // Viewport
    int viewport_width_;
    int viewport_height_;

    // Configuration
    Config config_;

    // Initial state (for reset)
    glm::vec3 initial_target_;
    float initial_distance_;
};

} // namespace visualization
