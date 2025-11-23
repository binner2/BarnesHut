#version 450 core

// Input vertex attributes
layout(location = 0) in vec3 in_position;
layout(location = 1) in vec4 in_color;

// Output to fragment shader
out vec4 frag_color;
out float frag_distance;

// Uniform variables
uniform mat4 u_view_projection;
uniform vec3 u_camera_position;
uniform float u_point_size;
uniform float u_point_scale;  // Size attenuation factor

void main() {
    // Transform position to clip space
    vec4 world_position = vec4(in_position, 1.0);
    gl_Position = u_view_projection * world_position;

    // Calculate distance from camera for size attenuation
    frag_distance = length(in_position - u_camera_position);

    // Apply size attenuation based on distance
    // Closer particles appear larger, distant ones smaller
    float attenuation = u_point_scale / (1.0 + 0.01 * frag_distance);
    gl_PointSize = u_point_size * attenuation;

    // Clamp point size to reasonable range
    gl_PointSize = clamp(gl_PointSize, 1.0, 64.0);

    // Pass color to fragment shader
    frag_color = in_color;
}
