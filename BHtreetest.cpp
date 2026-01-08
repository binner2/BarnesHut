/**
 * Modern Barnes-Hut N-Body Simulation
 * Modernized from 2004 legacy code to C++20
 *
 * Features:
 * - Modern C++20 features (concepts, ranges, span, etc.)
 * - Smart pointers for memory safety
 * - CPU parallelization with OpenMP
 * - RAII and exception safety
 * - std::optional for error handling
 */

#include "file.h"
#include "tree.h"
#include <iostream>
#include <string>
#include <vector>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace barnes_hut;

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <filename> <theta> <particles_per_leaf>\n"
              << "  filename: Input file with particle data\n"
              << "  theta: Barnes-Hut opening angle (e.g., 0.5)\n"
              << "  particles_per_leaf: Max particles in leaf node (e.g., 10)\n"
              << "\nExample:\n"
              << "  " << program_name << " data.dat 0.5 10\n";
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string filename = argv[1];
    const Real theta = std::stod(argv[2]);
    const Index particles_per_leaf = std::stoull(argv[3]);

    std::cout << "=== Modern Barnes-Hut N-Body Simulation ===\n\n";

    // Read configuration
    auto config_opt = read_config_file(filename);
    if (!config_opt) {
        std::cerr << "Error: Failed to read configuration from " << filename << "\n";
        return 2;
    }

    auto config = *config_opt;
    config.theta = theta;
    config.particles_per_leaf = particles_per_leaf;

    // Read particle data
    auto particles_opt = read_particle_file(filename, config);
    if (!particles_opt) {
        std::cerr << "Error: Failed to read particle data from " << filename << "\n";
        return 3;
    }

    auto particles = std::move(*particles_opt);

    std::cout << "Starting Barnes-Hut simulation for " << particles.size() << " particles\n"
              << "  Time: " << config.start_time << " -> " << config.end_time
              << " (dt=" << config.time_step << ")\n"
              << "  Theta: " << theta << "\n"
              << "  Max particles per leaf: " << particles_per_leaf << "\n";

#ifdef _OPENMP
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "  OpenMP enabled with " << omp_get_num_threads() << " threads\n";
    }
#else
    std::cout << "  OpenMP not enabled (serial execution)\n";
#endif

    std::cout << "\n";

    // Create Barnes-Hut tree
    BarnesHutTree tree(particles, config.time_step, theta, particles_per_leaf);

    // Simulation loop
    Index step = 0;
    Real current_time = config.start_time;

    const Real output_interval = (config.end_time - config.start_time) / 10.0;
    Real next_output_time = config.start_time + output_interval;

    Timer simulation_timer;

    while (current_time < config.end_time) {
        // Perform one simulation step
        tree.simulation_step();

        current_time += config.time_step;
        step++;

        // Print statistics
        const auto& stats = tree.get_statistics();

        if constexpr (ENABLE_TIMING) {
            std::cout << "Step " << std::setw(4) << step
                      << " | Time: " << std::fixed << std::setprecision(3) << current_time
                      << " | Load: " << std::setprecision(4) << stats.time_load << "s"
                      << " | Upward: " << stats.time_upward << "s"
                      << " | Force: " << stats.time_force << "s"
                      << " | Total: " << stats.time_total << "s"
                      << " | Direct: " << stats.direct_force_count
                      << " | P-C: " << stats.particle_cell_interactions
                      << "\n";
        }

        // Write snapshot if needed
        if (current_time >= next_output_time) {
            std::ostringstream message;
            message << tree.get_statistics_string();

            write_particle_forces(
                particles,
                message.str(),
                theta,
                particles_per_leaf
            );

            next_output_time += output_interval;
        }

        // Clear tree for next iteration
        tree.clear_tree();
    }

    const double total_simulation_time = simulation_timer.elapsed();

    std::cout << "\n=== Simulation Complete ===\n"
              << "Total steps: " << step << "\n"
              << "Total time: " << total_simulation_time << " seconds\n"
              << "Average time per step: " << (total_simulation_time / step) << " seconds\n";

    return 0;
}
