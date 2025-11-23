#include "tree.h"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <limits>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace barnes_hut {

BarnesHutTree::BarnesHutTree(std::span<Particle> particles, Real timestep, Real theta, Index max_particles_per_leaf)
    : particles_(particles)
    , dt_(timestep)
    , theta_(theta)
    , max_particles_per_leaf_(max_particles_per_leaf)
    , root_(std::make_unique<Node>())
    , current_node_index_(0)
    , max_tree_level_(0) {

    // Initialize particle IDs
    for (Index i = 0; i < particles_.size(); ++i) {
        particles_[i].set_id(i);
    }

    // Pre-allocate node pool (estimate: ~3x number of particles)
    node_pool_.reserve(particles_.size() * 3);
}

void BarnesHutTree::simulation_step() {
    Timer total_timer;
    stats_ = Statistics{};

    // Build tree
    Timer load_timer;
    build_tree();
    stats_.time_load = load_timer.elapsed();

    // Compute mass distribution
    Timer upward_timer;
    compute_mass_distribution();
    stats_.time_upward = upward_timer.elapsed();

    // Calculate forces
    Timer force_timer;
    #ifdef _OPENMP
    calculate_forces_parallel();
    #else
    calculate_forces();
    #endif
    stats_.time_force = force_timer.elapsed();

    // Integrate particles
    integrate_particles();

    stats_.time_total = total_timer.elapsed();
    stats_.nodes_used = current_node_index_;
    stats_.nodes_available = node_pool_.size();
}

void BarnesHutTree::clear_tree() {
    root_->reset();
    reset_node_pool();
}

void BarnesHutTree::find_bounding_box(Vector3D& center, Real& size) const {
    if (particles_.empty()) {
        center = Vector3D{0.0};
        size = 1.0;
        return;
    }

    Vector3D min_pos = particles_[0].position();
    Vector3D max_pos = particles_[0].position();

    // Find bounding box
    for (const auto& particle : particles_) {
        const auto& pos = particle.position();
        for (int dim = 0; dim < NDIM; ++dim) {
            min_pos[dim] = std::min(min_pos[dim], pos[dim]);
            max_pos[dim] = std::max(max_pos[dim], pos[dim]);
        }
    }

    // Calculate center and size
    for (int dim = 0; dim < NDIM; ++dim) {
        const Real extent = max_pos[dim] - min_pos[dim];
        center[dim] = min_pos[dim] + extent * 0.5;
        size = std::max(size, extent);
    }

    // Add small margin to ensure all particles are inside
    size = std::ceil(size) + 1.0;
}

void BarnesHutTree::build_tree() {
    // Find bounding box
    Vector3D center{0.0};
    Real size = 0.0;
    find_bounding_box(center, size);

    // Setup root node
    root_->reset();
    root_->geo_center = center;
    root_->size = size;
    root_->level = 0;
    root_->index = 0;

    // Insert all particles
    for (Index i = 0; i < particles_.size(); ++i) {
        Node* current = root_.get();
        bool inserted = false;

        while (!inserted) {
            const int child_idx = which_child(particles_[i].position(), current);
            if (child_idx < 0) {
                // Error in calculation
                continue;
            }

            auto& child = current->children[child_idx];

            if (!child || child->type == NodeType::Empty) {
                // Add new leaf
                add_leaf(i, particles_[i], current, child_idx);
                current->particle_count++;
                inserted = true;
            }
            else if (child->type == NodeType::Leaf) {
                // Check if we can add to existing leaf
                if (child->particle_count < max_particles_per_leaf_) {
                    child->particle_count++;
                    child->particle_list.push_back(&particles_[i]);
                    particles_[i].set_parent(child.get());
                    current->particle_count++;
                    inserted = true;
                }
                else {
                    // Convert leaf to internal node
                    convert_leaf_to_internal(current, child_idx);
                    current->particle_count++;
                    current = child.get();
                }
            }
            else if (child->type == NodeType::Internal) {
                current->particle_count++;
                current = child.get();
            }
        }
    }
}

int BarnesHutTree::which_child(const Vector3D& position, const Node* node) const noexcept {
    int child_number = 0;

    for (int k = 0; k < NDIM; ++k) {
        if (position[k] >= node->geo_center[k]) {
            child_number += (1 << k);
        }
    }

    return (child_number >= 0 && child_number < NSUB) ? child_number : -1;
}

void BarnesHutTree::add_leaf(Index particle_idx, Particle& particle, Node* node, int child_idx) {
    node->type = NodeType::Internal;

    // Allocate new leaf from pool
    auto* new_leaf = allocate_node();
    new_leaf->level = node->level + 1;
    new_leaf->particle_count = 1;
    new_leaf->particle_list.push_back(&particle);
    particle.set_parent(new_leaf);

    new_leaf->size = node->size / 2.0;
    new_leaf->type = NodeType::Leaf;

    // Calculate leaf center
    for (int k = 0; k < NDIM; ++k) {
        if ((child_idx >> k) & 1) {
            new_leaf->geo_center[k] = node->geo_center[k] + new_leaf->size / 2.0;
        }
        else {
            new_leaf->geo_center[k] = node->geo_center[k] - new_leaf->size / 2.0;
        }
    }

    // Transfer ownership to parent
    node->children[child_idx] = std::unique_ptr<Node>(new_leaf);
    new_leaf->parent = node;

    // Update max tree level
    max_tree_level_ = std::max(max_tree_level_, new_leaf->level);
}

void BarnesHutTree::convert_leaf_to_internal(Node* node, int child_idx) {
    auto* old_leaf = node->children[child_idx].get();

    // Save particle list
    auto temp_particle_list = std::move(old_leaf->particle_list);
    old_leaf->particle_count = 0;
    old_leaf->type = NodeType::Internal;

    // Re-insert particles
    for (auto* particle : temp_particle_list) {
        const int new_child_idx = which_child(particle->position(), old_leaf);
        insert_particle(particle->id(), *particle, old_leaf);
    }
}

void BarnesHutTree::insert_particle(Index particle_idx, Particle& particle, Node* node) {
    const int child_idx = which_child(particle.position(), node);
    auto& child = node->children[child_idx];

    if (!child || child->type == NodeType::Empty) {
        add_leaf(particle_idx, particle, node, child_idx);
        node->particle_count++;
        child->parent = node;
    }
    else if (child->type == NodeType::Leaf) {
        if (child->particle_count < max_particles_per_leaf_) {
            child->particle_count++;
            child->particle_list.push_back(&particle);
            particle.set_parent(child.get());
            node->particle_count++;
        }
        else {
            convert_leaf_to_internal(node, child_idx);
            node->particle_count++;
            insert_particle(particle_idx, particle, child.get());
        }
    }
    else if (child->type == NodeType::Internal) {
        node->particle_count++;
        insert_particle(particle_idx, particle, child.get());
    }
}

void BarnesHutTree::compute_mass_distribution() {
    compute_center_of_mass(root_.get());
}

void BarnesHutTree::compute_center_of_mass(Node* node) {
    if (!node || node->type == NodeType::Empty) {
        return;
    }

    if (node->type == NodeType::Leaf) {
        // Calculate center of mass for leaf
        Vector3D cms{0.0};
        Real total_mass = 0.0;

        for (const auto* particle : node->particle_list) {
            cms += particle->mass() * particle->position();
            total_mass += particle->mass();
        }

        if (total_mass > 0.0) {
            node->mass_center = cms / total_mass;
            node->mass = total_mass;
        }
    }
    else if (node->type == NodeType::Internal) {
        // Recursively calculate for children
        Vector3D cms{0.0};
        Real total_mass = 0.0;

        for (auto& child : node->children) {
            if (child && child->type != NodeType::Empty) {
                compute_center_of_mass(child.get());
                cms += child->mass * child->mass_center;
                total_mass += child->mass;
            }
        }

        if (total_mass > 0.0) {
            node->mass_center = cms / total_mass;
            node->mass = total_mass;
        }
    }
}

void BarnesHutTree::calculate_forces() {
    // Reset forces
    for (auto& particle : particles_) {
        particle.force() = Vector3D{0.0};
    }

    // Calculate forces for each particle
    for (auto& particle : particles_) {
        for (const auto& child : root_->children) {
            if (child && child->type != NodeType::Empty) {
                interact(particle, child.get());
            }
        }
    }
}

void BarnesHutTree::calculate_forces_parallel() {
    // Reset forces
    #pragma omp parallel for
    for (Index i = 0; i < particles_.size(); ++i) {
        particles_[i].force() = Vector3D{0.0};
    }

    // Calculate forces in parallel
    #pragma omp parallel for schedule(dynamic)
    for (Index i = 0; i < particles_.size(); ++i) {
        for (const auto& child : root_->children) {
            if (child && child->type != NodeType::Empty) {
                interact(particles_[i], child.get());
            }
        }
    }
}

bool BarnesHutTree::is_well_separated(const Particle& particle, const Node& node) const noexcept {
    const Real r_squared = particle.position().squared_distance(node.mass_center);
    const Real r = std::sqrt(r_squared + EPSILON_SQUARED);

    return (node.size / r) <= theta_;
}

void BarnesHutTree::interact(Particle& particle, const Node* node) {
    if (!node || node->type == NodeType::Empty) {
        return;
    }

    if (is_well_separated(particle, *node)) {
        // Use multipole approximation
        particle_cell_interaction(particle, *node);
    }
    else {
        // Need to go deeper
        if (node->type == NodeType::Internal) {
            for (const auto& child : node->children) {
                if (child && child->type != NodeType::Empty) {
                    interact(particle, child.get());
                }
            }
        }
        else if (node->type == NodeType::Leaf) {
            // Direct calculation with all particles in leaf
            for (auto* other_particle : node->particle_list) {
                direct_force_calculation(particle, *other_particle);
            }
        }
    }
}

void BarnesHutTree::particle_cell_interaction(Particle& particle, const Node& cell) {
    stats_.particle_cell_interactions++;

    const Real r_squared = particle.position().squared_distance(cell.mass_center);
    const Real r_cubed = (r_squared + EPSILON_SQUARED) * std::sqrt(r_squared + EPSILON_SQUARED);
    const Vector3D r_vec = particle.position() - cell.mass_center;

    particle.force() += -GRAVITY * particle.mass() * cell.mass / r_cubed * r_vec;
}

void BarnesHutTree::direct_force_calculation(Particle& p1, Particle& p2) {
    // Don't calculate force with itself
    if (p1.id() == p2.id()) {
        return;
    }

    stats_.direct_force_count++;

    const Vector3D r_vec = p1.position() - p2.position();
    const Real r_squared = r_vec.squared_magnitude();

    // Numerical stability with softening
    const Real r_cubed = (r_squared + EPSILON_SQUARED) * std::sqrt(r_squared + EPSILON_SQUARED);

    const Vector3D force = -GRAVITY * p1.mass() * p2.mass() / r_cubed * r_vec;
    p1.force() += force;
}

void BarnesHutTree::integrate_particles() {
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (Index i = 0; i < particles_.size(); ++i) {
        particles_[i].integrate(dt_);
    }
}

Node* BarnesHutTree::allocate_node() {
    if (current_node_index_ < node_pool_.size()) {
        auto* node = node_pool_[current_node_index_].get();
        node->reset();
        current_node_index_++;
        return node;
    }

    // Allocate new node
    auto new_node = std::make_unique<Node>();
    auto* node_ptr = new_node.get();
    node_pool_.push_back(std::move(new_node));
    current_node_index_++;
    return node_ptr;
}

void BarnesHutTree::reset_node_pool() noexcept {
    current_node_index_ = 0;
}

std::string BarnesHutTree::get_statistics_string() const {
    std::ostringstream oss;
    oss << "DirectForce: " << stats_.direct_force_count
        << "; ParticleCell: " << stats_.particle_cell_interactions
        << "; NodesUsed: " << stats_.nodes_used
        << "; NodesAvailable: " << stats_.nodes_available
        << "; TimeLoad: " << stats_.time_load
        << "; TimeUpward: " << stats_.time_upward
        << "; TimeForce: " << stats_.time_force
        << "; TimeTotal: " << stats_.time_total;
    return oss.str();
}

void BarnesHutTree::display_tree(const Node* node, std::ostream& os) const {
    if (!node) {
        node = root_.get();
    }

    if (!node) {
        os << "Tree is empty\n";
        return;
    }

    if (node->type == NodeType::Internal) {
        display_node(node, os);
        for (int i = 0; i < NSUB; ++i) {
            if (node->children[i] && node->children[i]->type != NodeType::Empty) {
                display_tree(node->children[i].get(), os);
            }
        }
    }
    else if (node->type == NodeType::Leaf) {
        display_node(node, os);
    }
}

void BarnesHutTree::display_node(const Node* node, std::ostream& os) const {
    os << std::fixed << std::setprecision(2);
    os << " Id=" << node->index
       << " L=" << node->level
       << " M=" << node->mass
       << " N=" << node->particle_count
       << " Geo=" << node->geo_center
       << " Size=" << node->size
       << " CMS=" << node->mass_center;

    switch (node->type) {
        case NodeType::Internal:
            os << " Type=Internal\n";
            break;
        case NodeType::Empty:
            os << " Type=Empty\n";
            break;
        case NodeType::Leaf:
            os << " Type=Leaf\n";
            for (Index i = 0; i < node->particle_list.size(); ++i) {
                os << "  Particle " << (i + 1) << " ID=" << node->particle_list[i]->id() << " ";
                node->particle_list[i]->display(os);
                os << "\n";
            }
            break;
    }
}

} // namespace barnes_hut
