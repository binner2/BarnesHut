#pragma once

#include "particle.h"
#include "vektor.h"
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

    // Main simulation step: kick-drift-kick leapfrog. The first call also
    // computes the initial forces.
    void simulation_step();

    // Rebuild the tree and compute forces at the current positions without
    // advancing time (used for validation and initial conditions).
    void compute_forces();

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

    // Release the current tree (nodes return to the pool). compute_forces()
    // does this itself; calling it between steps is harmless.
    void clear_tree();

    // Display tree (for debugging)
    void display_tree(const Node* node = nullptr, std::ostream& os = std::cout) const;

private:
    // Tree construction
    void find_bounding_box(Vector3D& center, Real& size) const;
    void build_tree();
    void insert_particle(Particle& particle, Node* node);
    [[nodiscard]] int which_child(const Vector3D& position, const Node* node) const noexcept;
    void add_leaf(Particle& particle, Node* node, int child_idx);
    void split_leaf(Node* leaf);

    // Tree traversal
    void compute_mass_distribution();
    void compute_center_of_mass(Node* node);

    // Force calculation
    // Interaction counters are accumulated per thread and reduced afterwards;
    // writing them straight into stats_ from the OpenMP loop is a data race
    // and serialises the threads on one cache line.
    struct InteractionCounts {
        Index direct = 0;
        Index particle_cell = 0;
    };

    void calculate_forces();
    void interact(Particle& particle, const Node* node, InteractionCounts& counts) const;
    [[nodiscard]] bool is_well_separated(const Particle& particle, const Node& node) const noexcept;
    void particle_cell_interaction(Particle& particle, const Node& cell) const;
    void direct_force_calculation(Particle& p1, const Particle& p2) const;

    // Integration
    void kick_particles(Real h);
    void drift_particles(Real dt);

    // Node management
    [[nodiscard]] Node* allocate_node();
    void reset_node_pool() noexcept;

    // Helper methods
    void display_node(const Node* node, std::ostream& os = std::cout) const;

    // Member variables
    std::span<Particle> particles_;
    Real dt_;
    Real theta_;
    Index max_particles_per_leaf_;

    std::unique_ptr<Node> root_;
    std::vector<std::unique_ptr<Node>> node_pool_;
    Index current_node_index_;

    Statistics stats_;
    Index max_tree_level_;
    bool forces_valid_ = false;
};

} // namespace barnes_hut
