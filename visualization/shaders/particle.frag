#version 450 core

// Input from vertex shader
in vec4 frag_color;
in float frag_distance;

// Output color
out vec4 out_color;

// Uniform variables
uniform float u_brightness;
uniform bool u_enable_glow;
uniform float u_alpha;

void main() {
    // Calculate circular point with smooth edges
    vec2 coord = gl_PointCoord * 2.0 - 1.0;
    float distance_from_center = length(coord);

    // Discard fragments outside the circle
    if (distance_from_center > 1.0) {
        discard;
    }

    // Create smooth circular particle with soft edges
    float alpha = 1.0 - smoothstep(0.5, 1.0, distance_from_center);

    // Apply glow effect if enabled
    vec3 final_color = frag_color.rgb * u_brightness;
    if (u_enable_glow) {
        // Add glow to center of particle
        float glow = exp(-distance_from_center * 3.0);
        final_color += vec3(glow * 0.3);
    }

    // Apply alpha and user-defined transparency
    alpha *= frag_color.a * u_alpha;

    // Fade distant particles slightly
    float distance_fade = 1.0 - smoothstep(50.0, 200.0, frag_distance);
    alpha *= mix(0.3, 1.0, distance_fade);

    // Output final color
    out_color = vec4(final_color, alpha);

    // Ensure we don't output pure black (looks better)
    if (out_color.a < 0.01) {
        discard;
    }
}
