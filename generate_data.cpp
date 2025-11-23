/**
 * Modern Data Generator for Barnes-Hut Simulation
 * Generates random particle data for testing
 */

#include "file.h"
#include <iostream>
#include <string>

using namespace barnes_hut;

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <filename> <num_particles> <t_start> <t_end> <dt>\n"
              << "  filename: Output file name\n"
              << "  num_particles: Number of particles to generate\n"
              << "  t_start: Simulation start time\n"
              << "  t_end: Simulation end time\n"
              << "  dt: Time step\n"
              << "\nExample:\n"
              << "  " << program_name << " test.dat 1000 0.0 1.0 0.01\n";
}

int main(int argc, char* argv[]) {
    std::cout << "=== Modern Barnes-Hut Data Generator ===\n\n";

    if (argc != 6) {
        print_usage(argv[0]);
        return 1;
    }

    const std::string filename = argv[1];
    const Index num_particles = std::stoull(argv[2]);
    const Real t_start = std::stod(argv[3]);
    const Real t_end = std::stod(argv[4]);
    const Real dt = std::stod(argv[5]);

    SimulationConfig config;
    config.particle_count = num_particles;
    config.start_time = t_start;
    config.end_time = t_end;
    config.time_step = dt;

    if (!config.is_valid()) {
        std::cerr << "Error: Invalid configuration parameters\n";
        std::cerr << "  Ensure: N > 0, t_end > t_start, dt > 0\n";
        return 2;
    }

    std::cout << "Generating data:\n"
              << "  Output file: " << filename << "\n"
              << "  Particles: " << num_particles << "\n"
              << "  Time range: " << t_start << " -> " << t_end << "\n"
              << "  Time step: " << dt << "\n\n";

    if (!generate_test_data(filename, config)) {
        std::cerr << "Error: Failed to generate data\n";
        return 3;
    }

    std::cout << "\nData generation successful!\n";
    return 0;
}
