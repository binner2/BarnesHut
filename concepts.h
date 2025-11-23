#pragma once

#include <concepts>
#include <type_traits>

namespace barnes_hut {

// Forward declarations
class Vector3D;
class Particle;
struct Node;

// ============================================================================
// Numeric Concepts
// ============================================================================

/// Concept for numeric types used in physics calculations
template<typename T>
concept Numeric = std::floating_point<T> || std::integral<T>;

/// Concept for scalar types (floating point only, for physics)
template<typename T>
concept Scalar = std::floating_point<T>;

/// Concept for index types
template<typename T>
concept IndexType = std::unsigned_integral<T> && sizeof(T) >= sizeof(std::size_t);

// ============================================================================
// Vector Concepts
// ============================================================================

/// Concept for 3D vector types
template<typename T>
concept Vector3DType = requires(T v, std::size_t i) {
    { v[i] } -> std::convertible_to<double>;
    { v.x() } -> std::convertible_to<double>;
    { v.y() } -> std::convertible_to<double>;
    { v.z() } -> std::convertible_to<double>;
    { v.squared_magnitude() } -> std::convertible_to<double>;
    { v + v } -> std::same_as<T>;
    { v - v } -> std::same_as<T>;
    { v * 1.0 } -> std::same_as<T>;
    { v / 1.0 } -> std::same_as<T>;
};

/// Concept for types that have a position
template<typename T>
concept HasPosition = requires(T obj) {
    { obj.position() } -> Vector3DType;
};

/// Concept for types that have a mass
template<typename T>
concept HasMass = requires(T obj) {
    { obj.mass() } -> std::convertible_to<double>;
};

// ============================================================================
// Physics Concepts
// ============================================================================

/// Concept for particle-like objects
template<typename T>
concept ParticleType = HasPosition<T> && HasMass<T> && requires(T p) {
    { p.velocity() } -> Vector3DType;
    { p.force() } -> Vector3DType;
    { p.id() } -> std::convertible_to<std::size_t>;
};

/// Concept for tree node types
template<typename T>
concept TreeNodeType = requires(T node) {
    { node.mass } -> std::convertible_to<double>;
    { node.mass_center } -> Vector3DType;
    { node.geo_center } -> Vector3DType;
    { node.size } -> std::convertible_to<double>;
    { node.particle_count } -> std::convertible_to<std::size_t>;
};

// ============================================================================
// Container Concepts
// ============================================================================

/// Concept for contiguous containers
template<typename T>
concept ContiguousContainer = requires(T container) {
    { container.data() };
    { container.size() } -> std::convertible_to<std::size_t>;
    { container.begin() };
    { container.end() };
};

/// Concept for particle containers
template<typename T>
concept ParticleContainer = ContiguousContainer<T> && requires(T container) {
    requires ParticleType<typename T::value_type>;
};

// ============================================================================
// Function Concepts
// ============================================================================

/// Concept for force calculation functions
template<typename F, typename P>
concept ForceCalculator = ParticleType<P> && requires(F func, P& p1, const P& p2) {
    { func(p1, p2) } -> std::same_as<void>;
};

/// Concept for integration schemes
template<typename F, typename P>
concept Integrator = ParticleType<P> && requires(F func, P& particle, double dt) {
    { func(particle, dt) } -> std::same_as<void>;
};

// ============================================================================
// Simulation Concepts
// ============================================================================

/// Concept for simulation time step
template<typename T>
concept TimeStepType = Scalar<T> && requires(T dt) {
    requires std::is_same_v<T, double> || std::is_same_v<T, float>;
};

/// Concept for types that can be simulated
template<typename T>
concept Simulatable = requires(T sim) {
    { sim.simulation_step() } -> std::same_as<void>;
};

// ============================================================================
// Memory Management Concepts
// ============================================================================

/// Concept for types that can be pooled
template<typename T>
concept Poolable = std::movable<T> && requires(T obj) {
    { obj.reset() } -> std::same_as<void>;
};

/// Concept for object pools
template<typename P, typename T>
concept ObjectPool = Poolable<T> && requires(P pool) {
    { pool.allocate() } -> std::convertible_to<T*>;
    { pool.reset() } -> std::same_as<void>;
};

// ============================================================================
// Utility Concepts
// ============================================================================

/// Concept for types that can be printed
template<typename T>
concept Printable = requires(T obj, std::ostream& os) {
    { os << obj } -> std::same_as<std::ostream&>;
};

/// Concept for types with statistics
template<typename T>
concept HasStatistics = requires(T obj) {
    { obj.get_statistics() };
    { obj.get_statistics_string() } -> std::convertible_to<std::string>;
};

// ============================================================================
// Compile-time checks
// ============================================================================

// Verify that our main types satisfy the concepts
static_assert(Vector3DType<Vector3D>, "Vector3D must satisfy Vector3DType");
static_assert(ParticleType<Particle>, "Particle must satisfy ParticleType");
static_assert(TreeNodeType<Node>, "Node must satisfy TreeNodeType");
static_assert(Poolable<Node>, "Node must be Poolable");

} // namespace barnes_hut
