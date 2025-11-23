#pragma once

#include "particle.h"
#include "vektor.h"
#include "object_pool.h"
#include <vector>
#include <memory>
#include <string>
#include <span>

namespace barnes_hut {

// Modern Barnes-Hut tree class with CPU parallelization support
class BarnesHutTree {
public:
    // Constructor
    BarnesHutTree(std::span<Particle> particles, Real timestep, Real theta, Index max_particles_per_leaf);

    // Rule of 5 - delete copy, default move
    ~BarnesHutTree() = default;
    BarnesHutTree(const BarnesHutTree&) = delete;
    BarnesHutTree& operator=(const BarnesHutTree&) = delete;
    BarnesHutTree(BarnesHutTree&&) noexcept = default;
    BarnesHutTree& operator=(BarnesHutTree&&) noexcept = default;

    // Main simulation step
    void simulation_step();

    // Get statistics
    struct Statistics {
        Index direct_force_count = 0;
        Index particle_cell_interactions = 0;
        Index nodes_used = 0;
        Index nodes_available = 0;
        double time_load = 0.0;
        double time_upward = 0.0;
        double time_force = 0.0;
        double time_total = 0.0;
    };

    [[nodiscard]] const Statistics& get_statistics() const noexcept { return stats_; }
    [[nodiscard]] std::string get_statistics_string() const;

    // Clear tree for next iteration
    void clear_tree();

    // Display tree (for debugging)
    void display_tree(const Node* node = nullptr, std::ostream& os = std::cout) const;

private:
    // Tree construction
    void find_bounding_box(Vector3D& center, Real& size) const;
    void build_tree();
    void insert_particle(Index particle_idx, Particle& particle, Node* node);
    [[nodiscard]] int which_child(const Vector3D& position, const Node* node) const noexcept;
    void add_leaf(Index particle_idx, Particle& particle, Node* node, int child_idx);
    void convert_leaf_to_internal(Node* node, int child_idx);

    // Tree traversal
    void compute_mass_distribution();
    void compute_center_of_mass(Node* node);

    // Force calculation
    void calculate_forces();
    void calculate_forces_parallel();  // Parallel version
    void interact(Particle& particle, const Node* node);
    [[nodiscard]] bool is_well_separated(const Particle& particle, const Node& node) const noexcept;
    void particle_cell_interaction(Particle& particle, const Node& cell);
    void direct_force_calculation(Particle& p1, Particle& p2);

    // Integration
    void integrate_particles();

    // Node management (using optimized object pool)
    [[nodiscard]] Node* allocate_node();
    void reset_node_pool() noexcept;

    // Pool efficiency statistics
    [[nodiscard]] double get_pool_efficiency() const noexcept;
    [[nodiscard]] std::size_t get_pool_memory_usage() const noexcept;

    // Helper methods
    void display_node(const Node* node, std::ostream& os = std::cout) const;

    // Member variables
    std::span<Particle> particles_;
    Real dt_;
    Real theta_;
    Index max_particles_per_leaf_;

    std::unique_ptr<Node> root_;
    ObjectPool<Node> node_pool_;  // Optimized object pool

    Statistics stats_;
    Index max_tree_level_;
};

} // namespace barnes_hut
