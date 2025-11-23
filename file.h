#pragma once

#include "particle.h"
#include "vektor.h"
#include <string>
#include <string_view>
#include <filesystem>
#include <optional>
#include <span>

namespace barnes_hut {

// Modern file I/O with exception safety and RAII
struct SimulationConfig {
    Index particle_count = 0;
    Real start_time = 0.0;
    Real end_time = 1.0;
    Real time_step = 0.01;
    Real theta = 0.5;
    Index particles_per_leaf = 1;

    [[nodiscard]] bool is_valid() const noexcept {
        return particle_count > 0 &&
               end_time > start_time &&
               time_step > 0.0 &&
               theta > 0.0 &&
               particles_per_leaf > 0;
    }
};

// Modern file reading with std::optional for error handling
[[nodiscard]] std::optional<SimulationConfig>
read_config_file(std::string_view filename);

// Read particle data from file
[[nodiscard]] std::optional<std::vector<Particle>>
read_particle_file(std::string_view filename, const SimulationConfig& config);

// Write particle positions to file
bool write_particle_positions(
    std::span<const Particle> particles,
    std::string_view message,
    Real theta,
    Index particles_per_leaf,
    std::string_view base_filename = "snapPOS"
);

// Write particle forces to file
bool write_particle_forces(
    std::span<const Particle> particles,
    std::string_view message,
    Real theta,
    Index particles_per_leaf,
    std::string_view base_filename = "snapFORCE"
);

// Generate random test data
bool generate_test_data(std::string_view filename, const SimulationConfig& config);

} // namespace barnes_hut
