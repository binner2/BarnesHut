#pragma once
/**
 * Point-cloud algorithms built on PointOctree.
 *
 * Three of them are "N-body shaped" kernel sums and use the Barnes-Hut
 * far-field approximation (traverse_bh):
 *   - kernel_density      : Gaussian kernel density estimate per point
 *   - mean_shift          : mode seeking / clustering
 *   - wlop_resample       : weighted locally optimal projection resampling
 * Setting theta = 0 (or exact = true) makes them evaluate the sums exactly
 * through the same traversal, which gives the reference for error studies.
 */
#include "point_cloud.h"
#include "point_octree.h"

#include <optional>
#include <vector>

namespace barnes_hut::pc {

// --- Classical octree-based preprocessing -----------------------------------

PointCloud voxel_downsample(const PointCloud& cloud, Real voxel_size);

struct NormalOptions {
    Index k = 20;                          // neighbours for PCA (used when radius <= 0)
    Real radius = 0.0;                     // radius neighbourhood instead of kNN
    std::optional<Vector3D> viewpoint;     // orient normals towards this point (camera origin)
};
// Fills cloud.normals; returns surface variation lambda0 / (l0+l1+l2) per point.
std::vector<Real> estimate_normals(PointCloud& cloud, const PointOctree& tree, const NormalOptions& opt = {});

// Indices of points whose mean kNN distance is within mean + std_ratio * stddev.
std::vector<Index> statistical_outlier_filter(const PointCloud& cloud, const PointOctree& tree,
                                              Index k = 20, Real std_ratio = 2.0);
std::vector<Index> radius_outlier_filter(const PointCloud& cloud, const PointOctree& tree,
                                         Real radius, Index min_neighbors);

// --- Barnes-Hut kernel sums --------------------------------------------------

struct BhCounters {
    Index far_evaluations = 0;    // node approximations used
    Index near_evaluations = 0;   // point-point evaluations
    double seconds = 0.0;
};

struct KernelDensityOptions {
    Real bandwidth = 0.02;        // h, metres
    Real theta = 0.5;             // opening angle; 0 = exact
    Real cutoff_sigmas = 4.0;     // ignore nodes farther than cutoff_sigmas * h (0 = none)
    Real kappa = 0.0;             // kernel-aware limit: approximate only cells of width <= kappa * h (0 = off)
};
struct KernelDensityResult {
    std::vector<Real> density;    // normalised Gaussian KDE, 1/m^3
    BhCounters counters;
};
KernelDensityResult kernel_density(const PointCloud& cloud, const PointOctree& tree,
                                   const KernelDensityOptions& opt);

struct MeanShiftOptions {
    Real bandwidth = 0.05;
    Real theta = 0.5;
    Real cutoff_sigmas = 4.0;
    Real kappa = 0.0;             // see KernelDensityOptions::kappa
    int max_iterations = 50;
    Real tolerance = 0.0;         // convergence step length; 0 -> 1e-3 * bandwidth
    Real merge_radius = 0.0;      // modes closer than this are one cluster; 0 -> bandwidth / 2
    Index min_cluster_size = 1;   // smaller clusters become label -1
};
struct MeanShiftResult {
    std::vector<int> labels;      // -1 = noise / too small
    std::vector<Vector3D> modes;  // one per label
    std::vector<Index> cluster_sizes;
    Index total_iterations = 0;
    BhCounters counters;
};
MeanShiftResult mean_shift(const PointCloud& cloud, const PointOctree& tree, const MeanShiftOptions& opt);

struct WlopOptions {
    Index target_count = 0;       // number of output points (0 -> 10 % of input)
    Real support_radius = 0.05;   // h of the WLOP kernel theta(r) = exp(-16 r^2 / h^2)
    Real repulsion_mu = 0.45;     // Huang et al. recommend 0 < mu < 0.5
    int iterations = 20;
    Real theta = 0.5;             // Barnes-Hut opening angle for both sums; 0 = exact
    Real kappa = 0.0;             // kernel-aware limit relative to the kernel sigma h/(4 sqrt 2) (0 = off)
    unsigned seed = 42;
    bool density_weights = true;  // the "W" in WLOP: compensate non-uniform input density
};
struct WlopResult {
    PointCloud cloud;
    std::vector<Real> mean_step_per_iteration;
    BhCounters counters;
};
WlopResult wlop_resample(const PointCloud& input, const WlopOptions& opt);

// --- Quality metrics ---------------------------------------------------------

struct NearestNeighborStats {
    Real mean = 0.0, stddev = 0.0, min = 0.0, max = 0.0;
    [[nodiscard]] Real coefficient_of_variation() const noexcept { return mean > 0.0 ? stddev / mean : 0.0; }
};
NearestNeighborStats nearest_neighbor_statistics(const PointCloud& cloud);

struct ErrorStats {
    Real max_relative = 0.0;      // max |a-b| / |b|
    Real mean_relative = 0.0;
    Real rms_relative = 0.0;
};
ErrorStats relative_error(std::span<const Real> approx, std::span<const Real> exact);

} // namespace barnes_hut::pc
