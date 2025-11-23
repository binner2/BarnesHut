#include "camera.h"
#include <glm/gtc/constants.hpp>
#include <algorithm>

namespace visualization {

Camera::Camera()
    : target_(0.0f)
    , position_(0.0f, 0.0f, 50.0f)
    , distance_(50.0f)
    , azimuth_(0.0f)
    , elevation_(glm::pi<float>() / 4.0f)  // 45 degrees
    , viewport_width_(800)
    , viewport_height_(600)
    , initial_target_(0.0f)
    , initial_distance_(50.0f)
{
    update_matrices();
}

void Camera::initialize(const glm::vec3& target, float distance) {
    target_ = target;
    distance_ = distance;
    initial_target_ = target;
    initial_distance_ = distance;

    azimuth_ = 0.0f;
    elevation_ = glm::pi<float>() / 4.0f;

    update();
}

void Camera::handle_mouse_move(float delta_x, float delta_y, bool left_button, bool right_button) {
    if (left_button) {
        // Orbit camera
        azimuth_ -= delta_x * config_.rotation_speed;
        elevation_ += delta_y * config_.rotation_speed;

        // Clamp elevation to prevent flipping
        const float epsilon = 0.01f;
        elevation_ = std::clamp(elevation_, epsilon, glm::pi<float>() - epsilon);

        update_position();
    } else if (right_button) {
        // Pan camera
        // Calculate pan direction in camera space
        glm::vec3 right = glm::normalize(glm::cross(position_ - target_, glm::vec3(0.0f, 1.0f, 0.0f)));
        glm::vec3 up = glm::normalize(glm::cross(right, position_ - target_));

        float pan_factor = distance_ * config_.pan_speed;
        target_ -= right * (delta_x * pan_factor);
        target_ += up * (delta_y * pan_factor);

        update_position();
    }
}

void Camera::handle_mouse_wheel(float delta) {
    // Zoom in/out by adjusting distance
    distance_ -= delta * config_.zoom_speed;
    distance_ = std::clamp(distance_, config_.min_distance, config_.max_distance);

    update_position();
}

void Camera::update() {
    update_position();
    update_matrices();
}

void Camera::set_viewport(int width, int height) {
    viewport_width_ = width;
    viewport_height_ = height;
    update_matrices();
}

void Camera::reset() {
    target_ = initial_target_;
    distance_ = initial_distance_;
    azimuth_ = 0.0f;
    elevation_ = glm::pi<float>() / 4.0f;
    update();
}

void Camera::update_position() {
    // Convert spherical coordinates to Cartesian
    float x = distance_ * std::sin(elevation_) * std::cos(azimuth_);
    float y = distance_ * std::cos(elevation_);
    float z = distance_ * std::sin(elevation_) * std::sin(azimuth_);

    position_ = target_ + glm::vec3(x, y, z);
    update_matrices();
}

void Camera::update_matrices() {
    // View matrix (look at target)
    view_matrix_ = glm::lookAt(
        position_,
        target_,
        glm::vec3(0.0f, 1.0f, 0.0f)  // Up vector
    );

    // Projection matrix (perspective)
    float aspect = static_cast<float>(viewport_width_) / static_cast<float>(viewport_height_);
    projection_matrix_ = glm::perspective(
        glm::radians(config_.fov),
        aspect,
        config_.near_plane,
        config_.far_plane
    );
}

} // namespace visualization
