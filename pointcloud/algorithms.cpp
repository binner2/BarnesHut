#include "algorithms.h"
#include "../stdinc.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <random>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace barnes_hut::pc {

namespace {

constexpr Real kPi = 3.14159265358979323846;

// Symmetric 3x3 eigen-decomposition by cyclic Jacobi. Eigenvalues ascending,
// eigenvectors as columns of V (V[k] is k-th eigenvector).
void eigen_symmetric_3x3(std::array<std::array<Real, 3>, 3> a, std::array<Real, 3>& w,
                         std::array<Vector3D, 3>& v) {
    std::array<std::array<Real, 3>, 3> V{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};
    for (int sweep = 0; sweep < 50; ++sweep) {
        const Real off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
        if (off < 1e-30) break;
        for (int p = 0; p < 3; ++p)
            for (int q = p + 1; q < 3; ++q) {
                if (std::abs(a[p][q]) < 1e-300) continue;
                const Real theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
                const Real t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
                const Real c = 1 / std::sqrt(t * t + 1), s = t * c;
                for (int k = 0; k < 3; ++k) {
                    const Real akp = a[k][p], akq = a[k][q];
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for (int k = 0; k < 3; ++k) {
                    const Real apk = a[p][k], aqk = a[q][k];
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for (int k = 0; k < 3; ++k) {
                    const Real vkp = V[k][p], vkq = V[k][q];
                    V[k][p] = c * vkp - s * vkq;
                    V[k][q] = s * vkp + c * vkq;
                }
            }
    }
    std::array<int, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](int i, int j) { return a[i][i] < a[j][j]; });
    for (int k = 0; k < 3; ++k) {
        w[k] = a[order[k]][order[k]];
        v[k] = Vector3D{V[0][order[k]], V[1][order[k]], V[2][order[k]]};
    }
}

struct VoxelKey {
    long long x, y, z;
    bool operator==(const VoxelKey& o) const noexcept { return x == o.x && y == o.y && z == o.z; }
};
struct VoxelHash {
    std::size_t operator()(const VoxelKey& k) const noexcept {
        return std::hash<long long>{}(k.x * 73856093LL ^ k.y * 19349663LL ^ k.z * 83492791LL);
    }
};

inline int thread_count() {
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

} // namespace

// ---------------------------------------------------------------------------

PointCloud voxel_downsample(const PointCloud& cloud, Real voxel_size) {
    struct Acc {
        Vector3D p{0.0}, n{0.0};
        std::array<Real, 3> c{0, 0, 0};
        Real s = 0.0;
        Index count = 0;
    };
    std::unordered_map<VoxelKey, Acc, VoxelHash> grid;
    grid.reserve(cloud.size() / 4 + 1);
    const Real inv = 1.0 / voxel_size;
    for (Index i = 0; i < cloud.size(); ++i) {
        const auto& p = cloud.points[i];
        const VoxelKey key{static_cast<long long>(std::floor(p.x() * inv)), static_cast<long long>(std::floor(p.y() * inv)),
                           static_cast<long long>(std::floor(p.z() * inv))};
        Acc& a = grid[key];
        a.p += p;
        if (cloud.has_normals()) a.n += cloud.normals[i];
        if (cloud.has_colors()) {
            a.c[0] += cloud.colors[i].r;
            a.c[1] += cloud.colors[i].g;
            a.c[2] += cloud.colors[i].b;
        }
        if (cloud.has_scalars()) a.s += cloud.scalars[i];
        ++a.count;
    }
    PointCloud out;
    out.reserve(grid.size());
    for (const auto& [key, a] : grid) {
        const Real inv_n = 1.0 / static_cast<Real>(a.count);
        out.points.push_back(a.p * inv_n);
        if (cloud.has_normals()) out.normals.push_back(a.n.normalized());
        if (cloud.has_colors())
            out.colors.push_back(Color{static_cast<std::uint8_t>(a.c[0] * inv_n), static_cast<std::uint8_t>(a.c[1] * inv_n),
                                       static_cast<std::uint8_t>(a.c[2] * inv_n)});
        if (cloud.has_scalars()) out.scalars.push_back(a.s * inv_n);
    }
    return out;
}

std::vector<Real> estimate_normals(PointCloud& cloud, const PointOctree& tree, const NormalOptions& opt) {
    const Index n = cloud.size();
    cloud.normals.assign(n, Vector3D{0.0});
    std::vector<Real> variation(n, 0.0);

#pragma omp parallel
    {
        std::vector<std::pair<Real, Index>> knn;
        std::vector<Index> nb;
#pragma omp for schedule(dynamic, 256)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii) {
            const Index i = static_cast<Index>(ii);
            nb.clear();
            if (opt.radius > 0.0) {
                tree.radius_search(cloud.points[i], opt.radius, nb);
            } else {
                tree.knn_search(cloud.points[i], opt.k, knn);
                for (const auto& [d2, j] : knn) nb.push_back(j);
            }
            if (nb.size() < 3) {
                cloud.normals[i] = Vector3D{0.0, 0.0, -1.0};
                continue;
            }
            Vector3D mean{0.0};
            for (const Index j : nb) mean += cloud.points[j];
            mean /= static_cast<Real>(nb.size());
            std::array<std::array<Real, 3>, 3> cov{};
            for (const Index j : nb) {
                const Vector3D d = cloud.points[j] - mean;
                for (int a = 0; a < 3; ++a)
                    for (int b = a; b < 3; ++b) cov[a][b] += d[a] * d[b];
            }
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < a; ++b) cov[a][b] = cov[b][a];
            std::array<Real, 3> w;
            std::array<Vector3D, 3> v;
            eigen_symmetric_3x3(cov, w, v);
            Vector3D normal = v[0].normalized();
            if (opt.viewpoint) {
                if (normal.dot(*opt.viewpoint - cloud.points[i]) < 0.0) normal = -normal;
            }
            cloud.normals[i] = normal;
            const Real sum = w[0] + w[1] + w[2];
            variation[i] = sum > 0.0 ? std::max<Real>(w[0], 0.0) / sum : 0.0;
        }
    }
    return variation;
}

std::vector<Index> statistical_outlier_filter(const PointCloud& cloud, const PointOctree& tree, Index k, Real std_ratio) {
    const Index n = cloud.size();
    std::vector<Real> mean_dist(n, 0.0);
#pragma omp parallel
    {
        std::vector<std::pair<Real, Index>> knn;
#pragma omp for schedule(dynamic, 256)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii) {
            const Index i = static_cast<Index>(ii);
            tree.knn_search(cloud.points[i], k + 1, knn);   // first hit is the point itself
            Real acc = 0.0;
            Index cnt = 0;
            for (const auto& [d2, j] : knn) {
                if (j == i) continue;
                acc += std::sqrt(d2);
                ++cnt;
            }
            mean_dist[i] = cnt ? acc / static_cast<Real>(cnt) : 0.0;
        }
    }
    const Real mu = std::accumulate(mean_dist.begin(), mean_dist.end(), 0.0) / static_cast<Real>(std::max<Index>(n, 1));
    Real var = 0.0;
    for (const Real d : mean_dist) var += (d - mu) * (d - mu);
    const Real sigma = std::sqrt(var / static_cast<Real>(std::max<Index>(n, 1)));
    const Real threshold = mu + std_ratio * sigma;
    std::vector<Index> keep;
    keep.reserve(n);
    for (Index i = 0; i < n; ++i)
        if (mean_dist[i] <= threshold) keep.push_back(i);
    return keep;
}

std::vector<Index> radius_outlier_filter(const PointCloud& cloud, const PointOctree& tree, Real radius, Index min_neighbors) {
    const Index n = cloud.size();
    std::vector<char> ok(n, 0);
#pragma omp parallel
    {
        std::vector<Index> nb;
#pragma omp for schedule(dynamic, 256)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii) {
            tree.radius_search(cloud.points[static_cast<Index>(ii)], radius, nb);
            ok[static_cast<Index>(ii)] = nb.size() > min_neighbors;   // nb contains the point itself
        }
    }
    std::vector<Index> keep;
    for (Index i = 0; i < n; ++i)
        if (ok[i]) keep.push_back(i);
    return keep;
}

// ---------------------------------------------------------------------------
// Barnes-Hut kernel density
// ---------------------------------------------------------------------------

KernelDensityResult kernel_density(const PointCloud& cloud, const PointOctree& tree, const KernelDensityOptions& opt) {
    Timer timer;
    const Index n = cloud.size();
    KernelDensityResult res;
    res.density.assign(n, 0.0);
    const Real h = opt.bandwidth;
    const Real inv_2h2 = 1.0 / (2.0 * h * h);
    const Real cutoff = opt.cutoff_sigmas > 0.0 ? opt.cutoff_sigmas * h : 0.0;
    const Real width_limit = opt.kappa > 0.0 ? opt.kappa * h : 0.0;
    const Real norm = 1.0 / (static_cast<Real>(n) * std::pow(2.0 * kPi, 1.5) * h * h * h);
    const auto points = tree.points();

    Index far_total = 0, near_total = 0;
#pragma omp parallel reduction(+ : far_total, near_total)
    {
#pragma omp for schedule(dynamic, 256)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii) {
            const Index i = static_cast<Index>(ii);
            const Vector3D& q = cloud.points[i];
            Real sum = 0.0;
            Index far = 0, near = 0;
            tree.traverse_bh(
                q, opt.theta, cutoff,
                [&](const OctreeNode& node) {
                    sum += node.weight * std::exp(-q.squared_distance(node.centroid) * inv_2h2);
                    ++far;
                },
                [&](Index j) {
                    sum += std::exp(-q.squared_distance(points[j]) * inv_2h2);
                    ++near;
                },
                width_limit);
            res.density[i] = sum * norm;
            far_total += far;
            near_total += near;
        }
    }
    res.counters.far_evaluations = far_total;
    res.counters.near_evaluations = near_total;
    res.counters.seconds = timer.elapsed();
    return res;
}

// ---------------------------------------------------------------------------
// Barnes-Hut mean shift
// ---------------------------------------------------------------------------

MeanShiftResult mean_shift(const PointCloud& cloud, const PointOctree& tree, const MeanShiftOptions& opt) {
    Timer timer;
    const Index n = cloud.size();
    MeanShiftResult res;
    res.labels.assign(n, -1);
    if (n == 0) return res;

    const Real h = opt.bandwidth;
    const Real inv_2h2 = 1.0 / (2.0 * h * h);
    const Real cutoff = opt.cutoff_sigmas > 0.0 ? opt.cutoff_sigmas * h : 0.0;
    const Real tol = opt.tolerance > 0.0 ? opt.tolerance : 1e-3 * h;
    const Real tol2 = tol * tol;
    const Real width_limit = opt.kappa > 0.0 ? opt.kappa * h : 0.0;
    const auto points = tree.points();

    std::vector<Vector3D> converged(n);
    Index far_total = 0, near_total = 0, iter_total = 0;

#pragma omp parallel reduction(+ : far_total, near_total, iter_total)
    {
#pragma omp for schedule(dynamic, 128)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii) {
            const Index i = static_cast<Index>(ii);
            Vector3D y = cloud.points[i];
            for (int it = 0; it < opt.max_iterations; ++it) {
                Vector3D num{0.0};
                Real den = 0.0;
                tree.traverse_bh(
                    y, opt.theta, cutoff,
                    [&](const OctreeNode& node) {
                        const Real w = node.weight * std::exp(-y.squared_distance(node.centroid) * inv_2h2);
                        num += node.centroid * w;
                        den += w;
                        ++far_total;
                    },
                    [&](Index j) {
                        const Real w = std::exp(-y.squared_distance(points[j]) * inv_2h2);
                        num += points[j] * w;
                        den += w;
                        ++near_total;
                    },
                    width_limit);
                ++iter_total;
                if (den <= 0.0) break;
                const Vector3D y_new = num / den;
                const Real step2 = y_new.squared_distance(y);
                y = y_new;
                if (step2 < tol2) break;
            }
            converged[i] = y;
        }
    }

    // Merge modes that lie within merge_radius of each other (union-find on a mode octree).
    const Real merge_r = opt.merge_radius > 0.0 ? opt.merge_radius : 0.5 * h;
    PointOctree mode_tree(converged, PointOctree::Params{32, 21, 0.0});
    std::vector<Index> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    auto find = [&](Index x) {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };
    std::vector<Index> nb;
    for (Index i = 0; i < n; ++i) {
        mode_tree.radius_search(converged[i], merge_r, nb);
        for (const Index j : nb) {
            const Index a = find(i), b = find(j);
            if (a != b) parent[std::max(a, b)] = std::min(a, b);
        }
    }
    std::unordered_map<Index, int> root_to_label;
    std::vector<Vector3D> mode_acc;
    std::vector<Index> sizes;
    for (Index i = 0; i < n; ++i) {
        const Index r = find(i);
        auto it = root_to_label.find(r);
        int label;
        if (it == root_to_label.end()) {
            label = static_cast<int>(mode_acc.size());
            root_to_label.emplace(r, label);
            mode_acc.push_back(Vector3D{0.0});
            sizes.push_back(0);
        } else {
            label = it->second;
        }
        res.labels[i] = label;
        mode_acc[label] += converged[i];
        ++sizes[label];
    }
    // Relabel by descending size, dropping small clusters.
    std::vector<int> order(mode_acc.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) { return sizes[a] > sizes[b]; });
    std::vector<int> remap(mode_acc.size(), -1);
    for (int rank = 0; rank < static_cast<int>(order.size()); ++rank) {
        const int old = order[rank];
        if (sizes[old] < opt.min_cluster_size) break;
        remap[old] = rank;
        res.modes.push_back(mode_acc[old] / static_cast<Real>(sizes[old]));
        res.cluster_sizes.push_back(sizes[old]);
    }
    for (Index i = 0; i < n; ++i) res.labels[i] = remap[res.labels[i]];

    res.total_iterations = iter_total;
    res.counters.far_evaluations = far_total;
    res.counters.near_evaluations = near_total;
    res.counters.seconds = timer.elapsed();
    return res;
}

// ---------------------------------------------------------------------------
// Barnes-Hut WLOP
// ---------------------------------------------------------------------------

WlopResult wlop_resample(const PointCloud& input, const WlopOptions& opt) {
    Timer timer;
    WlopResult res;
    const Index n = input.size();
    if (n == 0) return res;
    const Index m = opt.target_count > 0 ? std::min(opt.target_count, n) : std::max<Index>(1, n / 10);
    const Real h = opt.support_radius;
    const Real inv_h2_16 = 16.0 / (h * h);                  // theta(r) = exp(-16 r^2 / h^2), support ~ h
    const Real eps2 = 1e-12 * h * h;
    const Real width_limit = opt.kappa > 0.0 ? opt.kappa * h / (4.0 * std::sqrt(2.0)) : 0.0;

    // Input octree and WLOP density weights v_j = 1 + sum_k theta(|p_j - p_k|).
    PointOctree input_tree(input.points, PointOctree::Params{16, 21, 0.0});
    std::vector<Real> v(n, 1.0);
    if (opt.density_weights) {
#pragma omp parallel for schedule(dynamic, 256)
        for (std::ptrdiff_t jj = 0; jj < static_cast<std::ptrdiff_t>(n); ++jj) {
            const Index j = static_cast<Index>(jj);
            const Vector3D& p = input.points[j];
            Real acc = 0.0;
            input_tree.traverse_bh(
                p, opt.theta, h,
                [&](const OctreeNode& node) { acc += node.weight * std::exp(-p.squared_distance(node.centroid) * inv_h2_16); },
                [&](Index k) { acc += std::exp(-p.squared_distance(input.points[k]) * inv_h2_16); },
                width_limit);
            v[j] = acc;   // includes the self term (theta(0) = 1), so v_j >= 1
        }
    }
    // Attraction uses weights 1 / v_j; rebuild the input tree with those weights so that
    // far-field nodes carry the correct weighted centroid.
    std::vector<Real> inv_v(n);
    for (Index j = 0; j < n; ++j) inv_v[j] = 1.0 / v[j];
    input_tree.build(input.points, PointOctree::Params{16, 21, 0.0}, inv_v);

    // Initial sample: random subset of the input.
    std::mt19937_64 rng(opt.seed);
    std::vector<Index> perm(n);
    std::iota(perm.begin(), perm.end(), 0);
    std::shuffle(perm.begin(), perm.end(), rng);
    std::vector<Vector3D> x(m);
    for (Index i = 0; i < m; ++i) x[i] = input.points[perm[i]];

    std::vector<Vector3D> x_next(m);
    std::vector<Real> w(m, 1.0);
    Index far_total = 0, near_total = 0;

    for (int it = 0; it < opt.iterations; ++it) {
        PointOctree sample_tree(x, PointOctree::Params{16, 21, 0.0});

        // Repulsion density weights w_i = 1 + sum_{i' != i} theta(|x_i - x_i'|).
#pragma omp parallel for schedule(dynamic, 256)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(m); ++ii) {
            const Index i = static_cast<Index>(ii);
            Real acc = 0.0;
            sample_tree.traverse_bh(
                x[i], opt.theta, h,
                [&](const OctreeNode& node) { acc += node.weight * std::exp(-x[i].squared_distance(node.centroid) * inv_h2_16); },
                [&](Index k) { acc += std::exp(-x[i].squared_distance(x[k]) * inv_h2_16); },
                width_limit);
            w[i] = acc;
        }
        PointOctree sample_tree_w(x, PointOctree::Params{16, 21, 0.0}, w);

        Real step_sum = 0.0;
#pragma omp parallel for schedule(dynamic, 128) reduction(+ : step_sum, far_total, near_total)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(m); ++ii) {
            const Index i = static_cast<Index>(ii);
            const Vector3D& xi = x[i];

            // Attraction: sum_j p_j alpha_ij / sum_j alpha_ij, alpha_ij = theta(r) / (r v_j).
            Vector3D a_num{0.0};
            Real a_den = 0.0;
            input_tree.traverse_bh(
                xi, opt.theta, h,
                [&](const OctreeNode& node) {
                    const Real r = std::sqrt(xi.squared_distance(node.centroid) + eps2);
                    const Real a = node.weight * std::exp(-r * r * inv_h2_16) / r;
                    a_num += node.centroid * a;
                    a_den += a;
                    ++far_total;
                },
                [&](Index j) {
                    const Real r = std::sqrt(xi.squared_distance(input.points[j]) + eps2);
                    const Real a = inv_v[j] * std::exp(-r * r * inv_h2_16) / r;
                    a_num += input.points[j] * a;
                    a_den += a;
                    ++near_total;
                },
                width_limit);

            // Repulsion: mu * sum_{i'} (x_i - x_i') beta_ii' / sum beta_ii', beta = theta(r) w_i' / r.
            Vector3D r_num{0.0};
            Real r_den = 0.0;
            sample_tree_w.traverse_bh(
                xi, opt.theta, h,
                [&](const OctreeNode& node) {
                    const Real r = std::sqrt(xi.squared_distance(node.centroid) + eps2);
                    const Real b = node.weight * std::exp(-r * r * inv_h2_16) / r;
                    r_num += (xi - node.centroid) * b;
                    r_den += b;
                    ++far_total;
                },
                [&](Index k) {
                    if (k == i) return;
                    const Real r = std::sqrt(xi.squared_distance(x[k]) + eps2);
                    const Real b = w[k] * std::exp(-r * r * inv_h2_16) / r;
                    r_num += (xi - x[k]) * b;
                    r_den += b;
                    ++near_total;
                },
                width_limit);

            Vector3D target = a_den > 0.0 ? a_num / a_den : xi;
            if (r_den > 0.0) target += r_num * (opt.repulsion_mu / r_den);
            x_next[i] = target;
            step_sum += std::sqrt(target.squared_distance(xi));
        }
        x.swap(x_next);
        res.mean_step_per_iteration.push_back(step_sum / static_cast<Real>(m));
    }

    res.cloud.points = std::move(x);
    res.counters.far_evaluations = far_total;
    res.counters.near_evaluations = near_total;
    res.counters.seconds = timer.elapsed();
    return res;
}

// ---------------------------------------------------------------------------
// Metrics
// ---------------------------------------------------------------------------

NearestNeighborStats nearest_neighbor_statistics(const PointCloud& cloud) {
    NearestNeighborStats s;
    const Index n = cloud.size();
    if (n < 2) return s;
    PointOctree tree(cloud.points, PointOctree::Params{16, 21, 0.0});
    std::vector<Real> d(n, 0.0);
#pragma omp parallel
    {
        std::vector<std::pair<Real, Index>> knn;
#pragma omp for schedule(dynamic, 256)
        for (std::ptrdiff_t ii = 0; ii < static_cast<std::ptrdiff_t>(n); ++ii) {
            tree.knn_search(cloud.points[static_cast<Index>(ii)], 2, knn);
            d[static_cast<Index>(ii)] = knn.size() > 1 ? std::sqrt(knn[1].first) : 0.0;
        }
    }
    s.mean = std::accumulate(d.begin(), d.end(), 0.0) / static_cast<Real>(n);
    Real var = 0.0;
    for (const Real x : d) var += (x - s.mean) * (x - s.mean);
    s.stddev = std::sqrt(var / static_cast<Real>(n));
    const auto [mn, mx] = std::minmax_element(d.begin(), d.end());
    s.min = *mn;
    s.max = *mx;
    return s;
}

ErrorStats relative_error(std::span<const Real> approx, std::span<const Real> exact) {
    ErrorStats e;
    const Index n = std::min(approx.size(), exact.size());
    if (n == 0) return e;
    Real sum = 0.0, sum2 = 0.0;
    for (Index i = 0; i < n; ++i) {
        const Real denom = std::abs(exact[i]) > 0.0 ? std::abs(exact[i]) : 1.0;
        const Real r = std::abs(approx[i] - exact[i]) / denom;
        e.max_relative = std::max(e.max_relative, r);
        sum += r;
        sum2 += r * r;
    }
    e.mean_relative = sum / static_cast<Real>(n);
    e.rms_relative = std::sqrt(sum2 / static_cast<Real>(n));
    return e;
}

} // namespace barnes_hut::pc
