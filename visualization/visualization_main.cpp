/**
 * @file visualization_main.cpp
 * @brief Main entry point for Barnes-Hut simulation with CUDA-accelerated visualization
 *
 * This application combines:
 * - Barnes-Hut N-body simulation (CPU or GPU)
 * - CUDA-accelerated particle rendering
 * - OpenGL 4.5+ modern graphics pipeline
 * - Real-time interactive visualization
 */

#include <iostream>
#include <iomanip>
#include <GLFW/glfw3.h>

#include "file.h"
#include "tree.h"
#include "window_manager.h"
#include "camera.h"
#include "gl_renderer.h"
#include "cuda_renderer.cuh"

using namespace barnes_hut;
using namespace visualization;

/**
 * @brief Application state and configuration
 */
struct AppState {
    bool paused = false;
    bool single_step = false;
    bool show_help = false;
    int simulation_speed = 1;  // Steps per frame
    float time_scale = 1.0f;

    struct {
        bool color_by_velocity = true;
        bool color_by_force = false;
        bool enable_glow = true;
        float point_size = 3.0f;
    } render_config;
};

/**
 * @brief Print keyboard controls
 */
void print_controls() {
    std::cout << "\n=== Keyboard Controls ===\n"
              << "  Space     : Pause/Resume simulation\n"
              << "  S         : Single step (when paused)\n"
              << "  R         : Reset camera\n"
              << "  V         : Toggle velocity/force coloring\n"
              << "  G         : Toggle glow effect\n"
              << "  +/-       : Increase/Decrease particle size\n"
              << "  [/]       : Decrease/Increase simulation speed\n"
              << "  H         : Toggle help\n"
              << "  ESC       : Exit\n"
              << "\n=== Mouse Controls ===\n"
              << "  Left Drag : Rotate camera\n"
              << "  Right Drag: Pan camera\n"
              << "  Wheel     : Zoom in/out\n"
              << "=========================\n\n";
}

/**
 * @brief Display on-screen statistics
 */
void display_statistics(
    const BarnesHutTree::Statistics& sim_stats,
    const WindowManager& window,
    const AppState& app_state,
    Real current_time
) {
    if constexpr (ENABLE_TIMING) {
        std::cout << std::fixed << std::setprecision(3)
                  << "\rTime: " << std::setw(8) << current_time
                  << " | FPS: " << std::setw(6) << static_cast<int>(window.get_fps())
                  << " | Sim: " << std::setw(6) << sim_stats.time_total * 1000.0 << "ms"
                  << " | Paused: " << (app_state.paused ? "YES" : "NO ")
                  << " | Speed: " << app_state.simulation_speed << "x"
                  << std::flush;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        std::cout << "Usage: " << argv[0] << " <filename> <theta> <particles_per_leaf>\n"
                  << "  filename: Input file with particle data\n"
                  << "  theta: Barnes-Hut opening angle (e.g., 0.5)\n"
                  << "  particles_per_leaf: Max particles in leaf node (e.g., 10)\n"
                  << "\nExample:\n"
                  << "  " << argv[0] << " data.dat 0.5 10\n";
        return 1;
    }

    const std::string filename = argv[1];
    const Real theta = std::stod(argv[2]);
    const Index particles_per_leaf = std::stoull(argv[3]);

    std::cout << "=== Barnes-Hut N-Body Visualization ===\n\n";

    // ========================================================================
    // Load Simulation Data
    // ========================================================================

    std::cout << "Loading simulation data...\n";

    auto config_opt = read_config_file(filename);
    if (!config_opt) {
        std::cerr << "Error: Failed to read configuration from " << filename << "\n";
        return 2;
    }

    auto config = *config_opt;
    config.theta = theta;
    config.particles_per_leaf = particles_per_leaf;

    auto particles_opt = read_particle_file(filename, config);
    if (!particles_opt) {
        std::cerr << "Error: Failed to read particle data from " << filename << "\n";
        return 3;
    }

    auto particles = std::move(*particles_opt);
    const size_t particle_count = particles.size();

    std::cout << "Loaded " << particle_count << " particles\n";
    std::cout << "  Simulation time: " << config.start_time << " -> " << config.end_time << "\n";
    std::cout << "  Time step: " << config.time_step << "\n";
    std::cout << "  Theta: " << theta << "\n";
    std::cout << "  Particles per leaf: " << particles_per_leaf << "\n\n";

    // ========================================================================
    // Initialize Visualization
    // ========================================================================

    std::cout << "Initializing visualization system...\n";

    // Create window
    WindowManager window;
    WindowManager::Config window_config;
    window_config.title = "Barnes-Hut N-Body Simulation - CUDA Visualization";
    window_config.width = 1920;
    window_config.height = 1080;

    if (!window.initialize(window_config)) {
        std::cerr << "Failed to initialize window\n";
        return 4;
    }

    // Initialize OpenGL renderer
    GLRenderer gl_renderer;
    if (!gl_renderer.initialize(particle_count)) {
        std::cerr << "Failed to initialize OpenGL renderer\n";
        return 5;
    }

    // Initialize CUDA renderer
    CudaRenderer cuda_renderer;
    if (!cuda_renderer.initialize(particle_count)) {
        std::cerr << "Failed to initialize CUDA renderer\n";
        return 6;
    }

    // Register OpenGL VBO with CUDA for interop
    if (!cuda_renderer.register_gl_buffer(gl_renderer.get_vbo())) {
        std::cerr << "Failed to register OpenGL buffer with CUDA\n";
        return 7;
    }

    // Initialize camera
    Camera camera;
    int fb_width, fb_height;
    window.get_framebuffer_size(fb_width, fb_height);
    camera.set_viewport(fb_width, fb_height);

    // Calculate bounding box for initial camera position
    Vector3D min_pos(1e10), max_pos(-1e10);
    for (const auto& p : particles) {
        min_pos.x = std::min(min_pos.x, p.position().x);
        min_pos.y = std::min(min_pos.y, p.position().y);
        min_pos.z = std::min(min_pos.z, p.position().z);
        max_pos.x = std::max(max_pos.x, p.position().x);
        max_pos.y = std::max(max_pos.y, p.position().y);
        max_pos.z = std::max(max_pos.z, p.position().z);
    }
    Vector3D center = (min_pos + max_pos) * 0.5;
    Real extent = (max_pos - min_pos).length();
    camera.initialize(glm::vec3(center.x, center.y, center.z), extent * 1.5);

    std::cout << "Visualization initialized successfully\n";
    print_controls();

    // ========================================================================
    // Setup Input Callbacks
    // ========================================================================

    AppState app_state;

    WindowManager::InputCallbacks callbacks;

    callbacks.on_mouse_move = [&](float dx, float dy, bool left, bool right) {
        camera.handle_mouse_move(dx, dy, left, right);
    };

    callbacks.on_mouse_wheel = [&](float delta) {
        camera.handle_mouse_wheel(delta);
    };

    callbacks.on_key_press = [&](int key) {
        switch (key) {
            case GLFW_KEY_SPACE:
                app_state.paused = !app_state.paused;
                std::cout << "\nSimulation " << (app_state.paused ? "paused" : "resumed") << "\n";
                break;
            case GLFW_KEY_S:
                if (app_state.paused) {
                    app_state.single_step = true;
                }
                break;
            case GLFW_KEY_R:
                camera.reset();
                std::cout << "\nCamera reset\n";
                break;
            case GLFW_KEY_V:
                app_state.render_config.color_by_velocity = !app_state.render_config.color_by_velocity;
                app_state.render_config.color_by_force = !app_state.render_config.color_by_velocity;
                std::cout << "\nColoring by: " << (app_state.render_config.color_by_velocity ? "velocity" : "force") << "\n";
                break;
            case GLFW_KEY_G:
                app_state.render_config.enable_glow = !app_state.render_config.enable_glow;
                gl_renderer.get_config().enable_glow = app_state.render_config.enable_glow;
                std::cout << "\nGlow effect: " << (app_state.render_config.enable_glow ? "enabled" : "disabled") << "\n";
                break;
            case GLFW_KEY_EQUAL:  // + key
            case GLFW_KEY_KP_ADD:
                app_state.render_config.point_size += 0.5f;
                gl_renderer.get_config().point_size = app_state.render_config.point_size;
                std::cout << "\nPoint size: " << app_state.render_config.point_size << "\n";
                break;
            case GLFW_KEY_MINUS:
            case GLFW_KEY_KP_SUBTRACT:
                app_state.render_config.point_size = std::max(1.0f, app_state.render_config.point_size - 0.5f);
                gl_renderer.get_config().point_size = app_state.render_config.point_size;
                std::cout << "\nPoint size: " << app_state.render_config.point_size << "\n";
                break;
            case GLFW_KEY_LEFT_BRACKET:  // [
                app_state.simulation_speed = std::max(1, app_state.simulation_speed - 1);
                std::cout << "\nSimulation speed: " << app_state.simulation_speed << "x\n";
                break;
            case GLFW_KEY_RIGHT_BRACKET:  // ]
                app_state.simulation_speed = std::min(10, app_state.simulation_speed + 1);
                std::cout << "\nSimulation speed: " << app_state.simulation_speed << "x\n";
                break;
            case GLFW_KEY_H:
                app_state.show_help = !app_state.show_help;
                if (app_state.show_help) {
                    print_controls();
                }
                break;
        }
    };

    window.set_input_callbacks(callbacks);

    // ========================================================================
    // Create Barnes-Hut Tree
    // ========================================================================

    BarnesHutTree tree(particles, config.time_step, theta, particles_per_leaf);

    // ========================================================================
    // Main Simulation and Rendering Loop
    // ========================================================================

    std::cout << "Starting visualization loop...\n\n";

    Real current_time = config.start_time;
    Index step = 0;

    while (!window.should_close() && current_time < config.end_time) {
        window.begin_frame();

        // ====================================================================
        // Simulation Step(s)
        // ====================================================================

        if (!app_state.paused || app_state.single_step) {
            for (int i = 0; i < app_state.simulation_speed; ++i) {
                tree.simulation_step();
                current_time += config.time_step;
                step++;

                if (current_time >= config.end_time) {
                    break;
                }
            }

            app_state.single_step = false;
            tree.clear_tree();
        }

        // ====================================================================
        // Update Visualization
        // ====================================================================

        // Update camera
        camera.update();

        // Prepare CUDA renderer configuration
        CudaRenderer::Config cuda_config;
        cuda_config.color_by_velocity = app_state.render_config.color_by_velocity;
        cuda_config.color_by_force = app_state.render_config.color_by_force;
        cuda_config.point_size = app_state.render_config.point_size;

        // Update particle data on GPU via CUDA
        cuda_renderer.update_particles(particles.data(), particle_count, cuda_config);

        // Update OpenGL renderer camera
        gl_renderer.set_camera(
            camera.get_view_matrix(),
            camera.get_projection_matrix(),
            camera.get_position()
        );

        // Render particles
        gl_renderer.render(particle_count);

        // ====================================================================
        // Display Statistics
        // ====================================================================

        display_statistics(tree.get_statistics(), window, app_state, current_time);

        window.end_frame();
    }

    std::cout << "\n\n=== Simulation Complete ===\n";
    std::cout << "Total steps: " << step << "\n";
    std::cout << "Final time: " << current_time << "\n";

    // Cleanup
    cuda_renderer.cleanup();
    gl_renderer.cleanup();
    window.cleanup();

    return 0;
}
