#include "file.h"
#include "stdinc.h"
#include <fstream>
#include <sstream>
#include <iomanip>

namespace barnes_hut {

// Helper function for random number generation
namespace {
    inline double generate_random(double min_val, double max_val) {
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(min_val, max_val);
        return dist(gen);
    }
}  // anonymous namespace

std::optional<SimulationConfig> read_config_file(std::string_view filename) {
    std::ifstream infile(filename.data());
    if (!infile) {
        std::cerr << "Error: Could not open file: " << filename << "\n";
        return std::nullopt;
    }

    SimulationConfig config;

    // Read particle count
    if (!(infile >> config.particle_count) || config.particle_count <= 0) {
        std::cerr << "Error: Invalid particle count\n";
        return std::nullopt;
    }

    // Read start time
    if (!(infile >> config.start_time)) {
        std::cerr << "Error: Invalid start time\n";
        return std::nullopt;
    }

    // Read end time
    if (!(infile >> config.end_time) || config.end_time <= config.start_time) {
        std::cerr << "Error: Invalid end time (must be > start time)\n";
        return std::nullopt;
    }

    // Read time step
    if (!(infile >> config.time_step) || config.time_step <= 0.0) {
        std::cerr << "Error: Invalid time step\n";
        return std::nullopt;
    }

    return config;
}

std::optional<std::vector<Particle>> read_particle_file(
    std::string_view filename,
    const SimulationConfig& config) {

    std::ifstream infile(filename.data());
    if (!infile) {
        std::cerr << "Error: Could not open file: " << filename << "\n";
        return std::nullopt;
    }

    // Skip config data (already read)
    Index n;
    Real start_time, end_time, time_step;
    infile >> n >> start_time >> end_time >> time_step;

    if (n != config.particle_count) {
        std::cerr << "Error: Particle count mismatch\n";
        return std::nullopt;
    }

    std::vector<Particle> particles;
    particles.reserve(config.particle_count);

    for (Index i = 0; i < config.particle_count; ++i) {
        Real mass;
        Vector3D pos, vel;

        if (!(infile >> mass)) {
            std::cerr << "Error: Failed to read mass for particle " << i << "\n";
            return std::nullopt;
        }

        if (mass <= 0.0) {
            std::cerr << "Error: Invalid mass for particle " << i << " (must be > 0)\n";
            return std::nullopt;
        }

        // Read position
        for (int dim = 0; dim < NDIM; ++dim) {
            if (!(infile >> pos[dim])) {
                std::cerr << "Error: Failed to read position for particle " << i << "\n";
                return std::nullopt;
            }
        }

        // Read velocity
        for (int dim = 0; dim < NDIM; ++dim) {
            if (!(infile >> vel[dim])) {
                std::cerr << "Error: Failed to read velocity for particle " << i << "\n";
                return std::nullopt;
            }
        }

        particles.emplace_back(mass, pos, vel);
    }

    return particles;
}

bool write_particle_positions(
    std::span<const Particle> particles,
    std::string_view message,
    Real theta,
    Index particles_per_leaf,
    std::string_view base_filename) {

    static Index counter = 1;

    std::ostringstream filename;
    filename << base_filename << "_BH"
             << particles.size() / 1000 << "K"
             << "_theta" << std::fixed << std::setprecision(2) << theta
             << "_pLeaf" << particles_per_leaf
             << "_" << counter++ << ".dat";

    std::ofstream outfile(filename.str());
    if (!outfile) {
        std::cerr << "Error: Could not create output file: " << filename.str() << "\n";
        return false;
    }

    outfile << message << "\nNumber of particles = " << particles.size() << "\n\n";

    outfile << std::showpos << std::scientific << std::setprecision(6);
    for (const auto& particle : particles) {
        outfile << particle.position() << "\n";
    }

    std::cout << "Wrote positions to: " << filename.str() << "\n";
    return true;
}

bool write_particle_forces(
    std::span<const Particle> particles,
    std::string_view message,
    Real theta,
    Index particles_per_leaf,
    std::string_view base_filename) {

    static Index counter = 1;

    std::ostringstream filename;
    filename << base_filename << "_BH"
             << particles.size() / 1000 << "K"
             << "_theta" << std::fixed << std::setprecision(2) << theta
             << "_pLeaf" << particles_per_leaf
             << "_" << counter++ << ".dat";

    std::ofstream outfile(filename.str());
    if (!outfile) {
        std::cerr << "Error: Could not create output file: " << filename.str() << "\n";
        return false;
    }

    outfile << message << "\n" << particles.size() << "\n\n";

    outfile << std::showpos << std::scientific << std::setprecision(40);
    for (const auto& particle : particles) {
        outfile << particle.force() << "\n";
    }

    std::cout << "Wrote forces to: " << filename.str() << "\n";
    return true;
}

bool generate_test_data(std::string_view filename, const SimulationConfig& config) {
    if (!config.is_valid()) {
        std::cerr << "Error: Invalid configuration\n";
        return false;
    }

    std::ofstream outfile(filename.data());
    if (!outfile) {
        std::cerr << "Error: Could not create file: " << filename << "\n";
        return false;
    }

    // Write configuration
    outfile << config.particle_count << "\n"
            << config.start_time << "\n"
            << config.end_time << "\n"
            << config.time_step << "\n";

    // Generate random particles
    for (Index i = 0; i < config.particle_count; ++i) {
        // Mass between 5000 and 15000
        const Real mass = generate_random(5000.0, 15000.0);
        outfile << mass;

        // Position between 0 and 10
        for (int dim = 0; dim < NDIM; ++dim) {
            outfile << " " << generate_random(0.0, 10.0);
        }

        // Velocity between 0 and 100
        for (int dim = 0; dim < NDIM; ++dim) {
            outfile << " " << generate_random(0.0, 100.0);
        }

        outfile << "\n";
    }

    std::cout << "Generated test data in: " << filename << "\n";
    return true;
}

} // namespace barnes_hut
