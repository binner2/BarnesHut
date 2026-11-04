#include "point_octree.h"
#include "../stdinc.h"

#include <algorithm>
#include <array>
#include <numeric>
#include <queue>

namespace barnes_hut::pc {

PointOctree::PointOctree(std::span<const Vector3D> points, Params params, std::span<const Real> weights) {
    build(points, params, weights);
}

void PointOctree::build(std::span<const Vector3D> points, Params params, std::span<const Real> weights) {
    Timer timer;
    points_ = points;
    weights_ = weights;
    params_ = params;
    nodes_.clear();
    indices_.clear();
    stats_ = Statistics{};

    if (points_.empty()) return;
    if (!weights_.empty() && weights_.size() != points_.size())
        weights_ = {};

    // Bounding cube with a small margin so that max coordinates stay inside.
    Vector3D lo = points_[0], hi = points_[0];
    for (const auto& p : points_)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    Real extent = 0.0;
    Vector3D center{0.0};
    for (int k = 0; k < 3; ++k) {
        extent = std::max(extent, hi[k] - lo[k]);
        center[k] = 0.5 * (lo[k] + hi[k]);
    }
    Real half = 0.5 * extent * (1.0 + 1e-6) + 1e-9;
    if (half < 1e-9) half = 1e-6;

    indices_.resize(points_.size());
    std::iota(indices_.begin(), indices_.end(), 0u);
    nodes_.reserve(2 * points_.size() / std::max<Index>(1, params_.leaf_capacity) + 64);

    OctreeNode root;
    root.center = center;
    root.half = half;
    root.begin = 0;
    root.end = static_cast<std::uint32_t>(points_.size());
    root.level = 0;
    nodes_.push_back(root);

    std::vector<std::uint32_t> scratch(points_.size());
    build_recursive(0, 0, scratch);

    stats_.node_count = nodes_.size();
    stats_.leaf_count = static_cast<Index>(std::count_if(nodes_.begin(), nodes_.end(),
                                                         [](const OctreeNode& n) { return n.leaf; }));
    stats_.root_width = 2.0 * half;
    stats_.build_seconds = timer.elapsed();
}

void PointOctree::build_recursive(std::int32_t ni, int depth, std::vector<std::uint32_t>& scratch) {
    const std::uint32_t begin = nodes_[ni].begin;
    const std::uint32_t end = nodes_[ni].end;
    const std::uint32_t count = end - begin;
    const Vector3D center = nodes_[ni].center;
    const Real half = nodes_[ni].half;
    stats_.depth = std::max(stats_.depth, depth);

    const bool stop = count <= params_.leaf_capacity || depth >= params_.max_depth ||
                      (params_.min_half_size > 0.0 && 0.5 * half < params_.min_half_size);

    if (stop) {
        Vector3D acc{0.0};
        Real w_sum = 0.0;
        for (std::uint32_t i = begin; i < end; ++i) {
            const std::uint32_t idx = indices_[i];
            const Real w = weights_.empty() ? 1.0 : weights_[idx];
            acc += points_[idx] * w;
            w_sum += w;
        }
        OctreeNode& n = nodes_[ni];
        n.leaf = true;
        n.weight = w_sum;
        n.centroid = w_sum > 0.0 ? acc / w_sum : center;
        return;
    }

    // Counting sort of the range by octant (stable, O(count)).
    std::array<std::uint32_t, 8> counts{};
    for (std::uint32_t i = begin; i < end; ++i) {
        const Vector3D& p = points_[indices_[i]];
        int oct = 0;
        for (int k = 0; k < 3; ++k)
            if (p[k] >= center[k]) oct |= (1 << k);
        ++counts[oct];
    }
    std::array<std::uint32_t, 8> offsets{};
    std::uint32_t running = begin;
    for (int o = 0; o < 8; ++o) {
        offsets[o] = running;
        running += counts[o];
    }
    std::array<std::uint32_t, 8> cursor = offsets;
    for (std::uint32_t i = begin; i < end; ++i) {
        const std::uint32_t idx = indices_[i];
        const Vector3D& p = points_[idx];
        int oct = 0;
        for (int k = 0; k < 3; ++k)
            if (p[k] >= center[k]) oct |= (1 << k);
        scratch[cursor[oct]++] = idx;
    }
    std::copy(scratch.begin() + begin, scratch.begin() + end, indices_.begin() + begin);

    // Create children (all siblings before recursing keeps them contiguous).
    std::array<std::int32_t, 8> child_ids;
    child_ids.fill(-1);
    for (int o = 0; o < 8; ++o) {
        if (counts[o] == 0) continue;
        OctreeNode child;
        child.half = 0.5 * half;
        for (int k = 0; k < 3; ++k)
            child.center[k] = center[k] + (((o >> k) & 1) ? child.half : -child.half);
        child.begin = offsets[o];
        child.end = offsets[o] + counts[o];
        child.level = static_cast<std::uint8_t>(depth + 1);
        child_ids[o] = static_cast<std::int32_t>(nodes_.size());
        nodes_.push_back(child);
    }
    {
        OctreeNode& n = nodes_[ni];
        n.leaf = false;
        for (int o = 0; o < 8; ++o) n.children[o] = child_ids[o];
    }
    for (int o = 0; o < 8; ++o)
        if (child_ids[o] >= 0) build_recursive(child_ids[o], depth + 1, scratch);

    // Aggregate weighted centroid from children (bottom-up).
    Vector3D acc{0.0};
    Real w_sum = 0.0;
    for (int o = 0; o < 8; ++o) {
        if (child_ids[o] < 0) continue;
        const OctreeNode& c = nodes_[child_ids[o]];
        acc += c.centroid * c.weight;
        w_sum += c.weight;
    }
    OctreeNode& n = nodes_[ni];
    n.weight = w_sum;
    n.centroid = w_sum > 0.0 ? acc / w_sum : center;
}

void PointOctree::radius_search(const Vector3D& q, Real radius, std::vector<Index>& out) const {
    out.clear();
    if (nodes_.empty()) return;
    const Real r2 = radius * radius;
    std::vector<std::int32_t> stack;
    stack.reserve(64);
    stack.push_back(0);
    while (!stack.empty()) {
        const OctreeNode& n = nodes_[stack.back()];
        stack.pop_back();
        if (n.squared_distance_to_box(q) > r2) continue;
        if (n.leaf) {
            for (std::uint32_t i = n.begin; i < n.end; ++i) {
                const std::uint32_t idx = indices_[i];
                if (points_[idx].squared_distance(q) <= r2) out.push_back(idx);
            }
            continue;
        }
        for (const std::int32_t c : n.children)
            if (c >= 0) stack.push_back(c);
    }
}

void PointOctree::knn_search(const Vector3D& q, Index k, std::vector<std::pair<Real, Index>>& out) const {
    out.clear();
    if (nodes_.empty() || k == 0) return;

    using NodeEntry = std::pair<Real, std::int32_t>;             // (box distance^2, node)
    std::priority_queue<NodeEntry, std::vector<NodeEntry>, std::greater<NodeEntry>> frontier;
    std::priority_queue<std::pair<Real, Index>> best;             // max-heap on distance^2
    frontier.emplace(0.0, 0);

    while (!frontier.empty()) {
        const auto [d2_box, ni] = frontier.top();
        frontier.pop();
        if (best.size() >= k && d2_box >= best.top().first) break;
        const OctreeNode& n = nodes_[ni];
        if (n.leaf) {
            for (std::uint32_t i = n.begin; i < n.end; ++i) {
                const std::uint32_t idx = indices_[i];
                const Real d2 = points_[idx].squared_distance(q);
                if (best.size() < k) best.emplace(d2, idx);
                else if (d2 < best.top().first) {
                    best.pop();
                    best.emplace(d2, idx);
                }
            }
            continue;
        }
        for (const std::int32_t c : n.children)
            if (c >= 0) {
                const Real d2 = nodes_[c].squared_distance_to_box(q);
                if (best.size() < k || d2 < best.top().first) frontier.emplace(d2, c);
            }
    }
    out.reserve(best.size());
    while (!best.empty()) {
        out.push_back(best.top());
        best.pop();
    }
    std::reverse(out.begin(), out.end());
}

std::pair<Real, Index> PointOctree::nearest(const Vector3D& q) const {
    std::vector<std::pair<Real, Index>> tmp;
    knn_search(q, 1, tmp);
    if (tmp.empty()) return {std::numeric_limits<Real>::infinity(), npos};
    return tmp.front();
}

} // namespace barnes_hut::pc
