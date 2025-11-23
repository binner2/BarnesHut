#pragma once

#include "vektor.h"
#include <memory>
#include <cstdint>

namespace barnes_hut {

// Forward declaration
class TreeNode;

// Modern enum class (type-safe)
enum class NodeType : std::uint8_t {
    Empty = 0,
    Internal = 1,
    Leaf = 2
};

/**
 * @brief Cache-optimized tree node structure
 *
 * Layout optimized for cache performance:
 * - Frequently accessed fields (type, mass, centers) are placed first
 * - Hot path data fits within one cache line (64 bytes)
 * - Cold data (children pointers) placed at end
 */
struct alignas(64) Node {  // Cache-line aligned for better performance
    // ===== HOT DATA (frequently accessed during traversal) =====
    // First cache line (64 bytes)
    NodeType type = NodeType::Empty;        // 1 byte
    // 7 bytes padding here to align mass_center
    Vector3D mass_center{0.0};              // 24 bytes (aligned to 32 due to Vector3D alignment)
    Real mass = 0.0;                        // 8 bytes
    Real size = 0.0;                        // 8 bytes
    Vector3D geo_center{0.0};               // 24 bytes (aligned to 32)
    // Total hot data: ~72 bytes (spans into second cache line)

    // ===== WARM DATA (moderately accessed) =====
    Index particle_count = 0;               // 8 bytes
    Index level = 0;                        // 8 bytes
    Index index = 0;                        // 8 bytes
    Node* parent = nullptr;                 // 8 bytes (non-owning pointer)

    // ===== COLD DATA (less frequently accessed) =====
    std::vector<class Particle*> particle_list;  // For leaf nodes only
    std::array<std::unique_ptr<Node>, NSUB> children;  // Smart pointers!

    // Rule of 5 - default move, delete copy
    Node() = default;
    ~Node() = default;
    Node(const Node&) = delete;
    Node& operator=(const Node&) = delete;
    Node(Node&&) noexcept = default;
    Node& operator=(Node&&) noexcept = default;

    void reset() noexcept {
        index = 0;
        type = NodeType::Empty;
        geo_center = Vector3D{0.0};
        size = 0.0;
        mass_center = Vector3D{0.0};
        mass = 0.0;
        particle_list.clear();
        particle_count = 0;
        level = 0;
        for (auto& child : children) {
            child.reset();
        }
        parent = nullptr;
    }
};

// Modern Particle class with move semantics and const correctness
class Particle {
public:
    // Constructors
    Particle() noexcept
        : mass_(1.0)
        , position_(0.0)
        , velocity_(0.0)
        , force_(0.0)
        , id_(0)
        , parent_(nullptr) {}

    Particle(Real mass, Vector3D pos, Vector3D vel) noexcept
        : mass_(mass)
        , position_(pos)
        , velocity_(vel)
        , force_(0.0)
        , id_(0)
        , parent_(nullptr) {}

    // Rule of 5 - explicit defaults
    ~Particle() = default;
    Particle(const Particle&) = default;
    Particle& operator=(const Particle&) = default;
    Particle(Particle&&) noexcept = default;
    Particle& operator=(Particle&&) noexcept = default;

    // Const-correct getters (return const ref for efficiency)
    [[nodiscard]] Real mass() const noexcept { return mass_; }
    [[nodiscard]] const Vector3D& position() const noexcept { return position_; }
    [[nodiscard]] const Vector3D& velocity() const noexcept { return velocity_; }
    [[nodiscard]] const Vector3D& force() const noexcept { return force_; }
    [[nodiscard]] Index id() const noexcept { return id_; }
    [[nodiscard]] const Node* parent() const noexcept { return parent_; }

    // Non-const accessors for modification
    [[nodiscard]] Vector3D& position() noexcept { return position_; }
    [[nodiscard]] Vector3D& velocity() noexcept { return velocity_; }
    [[nodiscard]] Vector3D& force() noexcept { return force_; }

    // Setters
    void set_position(const Vector3D& pos) noexcept { position_ = pos; }
    void set_velocity(const Vector3D& vel) noexcept { velocity_ = vel; }
    void set_force(const Vector3D& f) noexcept { force_ = f; }
    void set_id(Index identity) noexcept { id_ = identity; }
    void set_parent(Node* p) noexcept { parent_ = p; }

    // Leapfrog integration for velocity and position
    void integrate(Real dt) noexcept {
        const Vector3D acceleration = force_ / mass_;

        // Leapfrog integration (kick-drift-kick)
        velocity_ += acceleration * (0.5 * dt);
        position_ += velocity_ * dt;
        velocity_ += acceleration * (0.5 * dt);
    }

    // Display particle information
    void display(std::ostream& os = std::cout) const {
        if (parent_) {
            os << "parent=" << parent_->index;
        }
        os << " mass=" << mass_ << " force=";
        force_.print(os);
        os << " pos=";
        position_.print(os);
    }

private:
    Real mass_;
    Vector3D position_;
    Vector3D velocity_;
    Vector3D force_;
    Index id_;
    Node* parent_;  // Non-owning pointer

    // Friend declaration for tree access
    friend class BarnesHutTree;
};

// Type alias for backward compatibility
using particlePtr = Particle*;

} // namespace barnes_hut
