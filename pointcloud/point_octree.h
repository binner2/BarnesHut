#pragma once
/**
 * Index-based (linear) octree over a static point set.
 *
 * Design differences from the simulation tree in ../tree.h:
 *  - nodes live in one contiguous std::vector, children are int32 indices;
 *  - every node owns a contiguous range [begin, end) of a permuted index
 *    array, so both leaves and internal nodes can enumerate their points;
 *  - no per-node heap allocations; the layout is GPU/SoA friendly.
 *
 * The octant convention (bit k set when p[k] >= center[k]) and the
 * Barnes-Hut opening criterion (cell width / distance < theta) follow the
 * original code so that results are directly comparable.
 */
#include "../vektor.h"
#include <cstdint>
#include <limits>
#include <span>
#include <utility>
#include <vector>

namespace barnes_hut::pc {

struct OctreeNode {
    Vector3D center{0.0};
    Real half = 0.0;                 // half edge length
    Vector3D centroid{0.0};          // weighted mean of points in the subtree
    Real weight = 0.0;               // sum of weights (= count for unit weights)
    std::uint32_t begin = 0;         // range in PointOctree::indices()
    std::uint32_t end = 0;
    std::int32_t children[8] = {-1, -1, -1, -1, -1, -1, -1, -1};
    std::uint8_t level = 0;
    bool leaf = true;

    [[nodiscard]] std::uint32_t count() const noexcept { return end - begin; }
    [[nodiscard]] Real width() const noexcept { return 2.0 * half; }

    // Squared distance from q to the node's cube (0 if inside).
    [[nodiscard]] Real squared_distance_to_box(const Vector3D& q) const noexcept {
        Real d2 = 0.0;
        for (int k = 0; k < 3; ++k) {
            const Real lo = center[k] - half, hi = center[k] + half;
            if (q[k] < lo) { const Real d = lo - q[k]; d2 += d * d; }
            else if (q[k] > hi) { const Real d = q[k] - hi; d2 += d * d; }
        }
        return d2;
    }
    [[nodiscard]] bool contains(const Vector3D& q) const noexcept {
        for (int k = 0; k < 3; ++k)
            if (q[k] < center[k] - half || q[k] >= center[k] + half) return false;
        return true;
    }
};

// Build parameters (kept outside the class so that it can serve as a default
// argument; GCC rejects nested aggregates with NSDMIs in that position).
struct OctreeParams {
    Index leaf_capacity = 16;    // split a node when it holds more points
    int max_depth = 21;          // guards against coincident points
    Real min_half_size = 0.0;    // stop splitting below this half size
};

class PointOctree {
public:
    using Params = OctreeParams;

    struct Statistics {
        Index node_count = 0;
        Index leaf_count = 0;
        int depth = 0;
        Real root_width = 0.0;
        double build_seconds = 0.0;
    };

    PointOctree() = default;
    PointOctree(std::span<const Vector3D> points, Params params = {}, std::span<const Real> weights = {});

    void build(std::span<const Vector3D> points, Params params = {}, std::span<const Real> weights = {});

    [[nodiscard]] const std::vector<OctreeNode>& nodes() const noexcept { return nodes_; }
    [[nodiscard]] const std::vector<std::uint32_t>& indices() const noexcept { return indices_; }
    [[nodiscard]] std::span<const Vector3D> points() const noexcept { return points_; }
    [[nodiscard]] const OctreeNode& root() const noexcept { return nodes_[0]; }
    [[nodiscard]] const Statistics& statistics() const noexcept { return stats_; }
    [[nodiscard]] bool empty() const noexcept { return nodes_.empty(); }

    // All points with |p - q| <= radius. Output is unsorted.
    void radius_search(const Vector3D& q, Real radius, std::vector<Index>& out) const;
    // k nearest points, sorted by ascending distance; pairs are (squared distance, index).
    void knn_search(const Vector3D& q, Index k, std::vector<std::pair<Real, Index>>& out) const;
    // Nearest point (squared distance, index); index = npos if empty.
    [[nodiscard]] std::pair<Real, Index> nearest(const Vector3D& q) const;

    /**
     * Barnes-Hut traversal. For every node the opening criterion
     *      width / |q - centroid| < theta
     * is tested; a node that passes is handed to far(node) once, otherwise
     * it is opened. Leaves that cannot be approximated hand each point index
     * to near(index). The node containing q is never approximated, so the
     * self term is always evaluated exactly. Nodes whose cube is farther than
     * cutoff from q are skipped entirely (cutoff <= 0 disables skipping).
     * theta <= 0 degenerates to an exact (brute-force through leaves) sum.
     *
     * far_width_limit > 0 adds a kernel-aware condition: a node is only
     * approximated when its width is also <= far_width_limit (e.g. kappa * h
     * for a Gaussian kernel of bandwidth h). The purely geometric criterion
     * bounds the angular size of a cell, which is the right quantity for a
     * 1/r potential but not for a kernel with an intrinsic length scale.
     */
    template <class FarFn, class NearFn>
    void traverse_bh(const Vector3D& q, Real theta, Real cutoff, FarFn&& far, NearFn&& near,
                     Real far_width_limit = 0.0) const {
        if (nodes_.empty()) return;
        const Real cutoff2 = cutoff > 0.0 ? cutoff * cutoff : std::numeric_limits<Real>::infinity();
        const Real width_limit = far_width_limit > 0.0 ? far_width_limit : std::numeric_limits<Real>::infinity();
        traverse_bh_impl(0, q, theta, cutoff2, width_limit, far, near);
    }

    static constexpr Index npos = std::numeric_limits<Index>::max();

private:
    void build_recursive(std::int32_t node_index, int depth, std::vector<std::uint32_t>& scratch);

    template <class FarFn, class NearFn>
    void traverse_bh_impl(std::int32_t ni, const Vector3D& q, Real theta, Real cutoff2, Real width_limit,
                          FarFn& far, NearFn& near) const {
        const OctreeNode& n = nodes_[ni];
        if (n.weight <= 0.0) return;
        if (n.squared_distance_to_box(q) > cutoff2) return;

        if (theta > 0.0 && !n.contains(q)) {
            const Real d2 = q.squared_distance(n.centroid);
            const Real w = n.width();
            if (w <= width_limit && w * w < theta * theta * d2) {
                far(n);
                return;
            }
        }
        if (n.leaf) {
            for (std::uint32_t i = n.begin; i < n.end; ++i) near(static_cast<Index>(indices_[i]));
            return;
        }
        for (const std::int32_t c : n.children)
            if (c >= 0) traverse_bh_impl(c, q, theta, cutoff2, width_limit, far, near);
    }

    std::span<const Vector3D> points_;
    std::span<const Real> weights_;
    Params params_;
    std::vector<OctreeNode> nodes_;
    std::vector<std::uint32_t> indices_;
    Statistics stats_;
};

} // namespace barnes_hut::pc
