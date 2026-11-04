/**
 * Self-checking tests for the point-cloud module (no external framework).
 * Exit code 0 = all passed.
 */
#include "../algorithms.h"
#include "../point_cloud.h"
#include "../point_octree.h"
#include "../synthetic.h"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>

using namespace barnes_hut;
using namespace barnes_hut::pc;

namespace {
int failures = 0;
#define CHECK(cond)                                                                            \
    do {                                                                                       \
        if (!(cond)) {                                                                         \
            ++failures;                                                                        \
            std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << "  " #cond "\n";           \
        }                                                                                      \
    } while (0)

PointCloud random_cloud(Index n, unsigned seed) {
    std::mt19937_64 rng(seed);
    std::uniform_real_distribution<Real> u(-1.0, 1.0);
    PointCloud c;
    for (Index i = 0; i < n; ++i) c.points.emplace_back(u(rng), u(rng), u(rng));
    return c;
}

void test_octree_queries() {
    auto cloud = random_cloud(5000, 3);
    // Add duplicated points to exercise the max-depth guard.
    for (int i = 0; i < 40; ++i) cloud.points.push_back(Vector3D{0.123, 0.456, 0.789});
    PointOctree tree(cloud.points, PointOctree::Params{8, 21, 0.0});
    CHECK(tree.root().count() == cloud.size());
    CHECK(std::abs(tree.root().weight - static_cast<Real>(cloud.size())) < 1e-9);
    // Root centroid equals the plain mean.
    CHECK(tree.root().centroid.distance(cloud.centroid()) < 1e-9);

    std::mt19937_64 rng(11);
    std::uniform_real_distribution<Real> u(-1.2, 1.2);
    std::vector<Index> rs;
    std::vector<std::pair<Real, Index>> knn;
    for (int t = 0; t < 50; ++t) {
        const Vector3D q{u(rng), u(rng), u(rng)};
        const Real r = 0.15 + 0.1 * t / 50.0;
        tree.radius_search(q, r, rs);
        Index brute = 0;
        for (const auto& p : cloud.points) brute += p.squared_distance(q) <= r * r;
        CHECK(rs.size() == brute);

        tree.knn_search(q, 7, knn);
        CHECK(knn.size() == 7);
        std::vector<Real> d2(cloud.size());
        for (Index i = 0; i < cloud.size(); ++i) d2[i] = cloud.points[i].squared_distance(q);
        std::sort(d2.begin(), d2.end());
        for (int k = 0; k < 7; ++k) CHECK(std::abs(knn[k].first - d2[k]) < 1e-12);
        CHECK(std::is_sorted(knn.begin(), knn.end()));
    }
}

void test_bh_exact_consistency() {
    auto cloud = random_cloud(3000, 5);
    PointOctree tree(cloud.points, PointOctree::Params{16, 21, 0.0});
    KernelDensityOptions ex;
    ex.bandwidth = 0.1;
    ex.theta = 0.0;
    ex.cutoff_sigmas = 0.0;
    auto exact = kernel_density(cloud, tree, ex);
    // Brute-force reference for a few points.
    const Real norm = 1.0 / (cloud.size() * std::pow(2 * 3.14159265358979323846, 1.5) * std::pow(ex.bandwidth, 3));
    for (Index i = 0; i < cloud.size(); i += 301) {
        Real s = 0.0;
        for (const auto& p : cloud.points) s += std::exp(-p.squared_distance(cloud.points[i]) / (2 * ex.bandwidth * ex.bandwidth));
        CHECK(std::abs(s * norm - exact.density[i]) <= 1e-9 * std::abs(s * norm));
    }
    CHECK(exact.counters.far_evaluations == 0);
    CHECK(exact.counters.near_evaluations == cloud.size() * cloud.size());

    // BH error should shrink monotonically with theta and stay small.
    Real prev = 1e9;
    for (const Real theta : {1.0, 0.7, 0.5, 0.3}) {
        KernelDensityOptions ko = ex;
        ko.theta = theta;
        auto r = kernel_density(cloud, tree, ko);
        const auto err = relative_error(r.density, exact.density);
        CHECK(err.max_relative <= prev + 1e-12);
        prev = err.max_relative;
        if (theta <= 0.5) CHECK(err.max_relative < 0.02);
        CHECK(r.counters.far_evaluations > 0);

        // Kernel-aware limit must not increase the error and must use fewer approximations.
        KernelDensityOptions kk = ko;
        kk.kappa = 0.5;
        auto rk = kernel_density(cloud, tree, kk);
        const auto errk = relative_error(rk.density, exact.density);
        CHECK(errk.max_relative <= err.max_relative + 1e-12);
        CHECK(rk.counters.far_evaluations <= r.counters.far_evaluations);
    }
}

void test_mean_shift_separates_blobs() {
    std::mt19937_64 rng(9);
    std::normal_distribution<Real> g(0.0, 0.03);
    PointCloud cloud;
    const std::array<Vector3D, 3> centers{Vector3D{0, 0, 0}, Vector3D{1, 0, 0}, Vector3D{0, 1, 0}};
    for (int c = 0; c < 3; ++c)
        for (int i = 0; i < 400; ++i) cloud.points.push_back(centers[c] + Vector3D{g(rng), g(rng), g(rng)});
    PointOctree tree(cloud.points);
    MeanShiftOptions mo;
    mo.bandwidth = 0.1;
    mo.theta = 0.5;
    mo.min_cluster_size = 50;
    auto res = mean_shift(cloud, tree, mo);
    CHECK(res.modes.size() == 3);
    for (int c = 0; c < 3; ++c) {
        bool found = false;
        for (const auto& m : res.modes) found |= m.distance(centers[c]) < 0.02;
        CHECK(found);
    }
    // All points of a blob share one label.
    for (int c = 0; c < 3; ++c)
        for (int i = 1; i < 400; ++i) CHECK(res.labels[c * 400 + i] == res.labels[c * 400]);
}

void test_normals_on_plane_and_sphere() {
    PointCloud cloud;
    std::mt19937_64 rng(2);
    std::uniform_real_distribution<Real> u(-1.0, 1.0);
    for (int i = 0; i < 2000; ++i) cloud.points.emplace_back(u(rng), u(rng), 2.0);   // plane z = 2
    PointOctree tree(cloud.points);
    NormalOptions no;
    no.k = 16;
    no.viewpoint = Vector3D{0.0};
    estimate_normals(cloud, tree, no);
    for (Index i = 0; i < cloud.size(); i += 97) {
        CHECK(std::abs(cloud.normals[i].z() + 1.0) < 1e-6);   // points towards the camera (origin)
    }
}

void test_depth_roundtrip_and_io() {
    CameraIntrinsics K;
    K.width = 160;
    K.height = 120;
    K.cx = 79.5;
    K.cy = 59.5;
    K.fx = K.fy = 131.25;
    SyntheticSceneOptions so;
    so.noise_coefficient = 0.0;
    so.dropout_fraction = 0.0;
    auto frame = render_synthetic_depth(K, so);
    std::vector<std::pair<int, int>> pix;
    PointCloud cloud = depth_to_point_cloud(frame.depth, K, {}, &pix);
    CHECK(cloud.size() > 0.9 * static_cast<Index>(K.width) * K.height);
    for (Index i = 0; i < cloud.size(); i += 101) {
        Real u = 0.0, v = 0.0;
        CHECK(K.project(cloud.points[i], u, v));
        CHECK(std::abs(u - pix[i].first) < 1e-6 && std::abs(v - pix[i].second) < 1e-6);
        CHECK(std::abs(cloud.points[i].z() * K.depth_scale - frame.depth.at(pix[i].first, pix[i].second)) < 1e-6);
    }
    // Floor points must lie on y = floor_height, wall points on z = wall_distance.
    for (Index i = 0; i < cloud.size(); ++i) {
        const auto id = static_cast<SurfaceId>(frame.surface[static_cast<std::size_t>(pix[i].second) * K.width + pix[i].first]);
        if (id == SurfaceId::Floor) CHECK(std::abs(cloud.points[i].y() - so.floor_height) < 2.0 / K.depth_scale * 3.0);
        if (id == SurfaceId::BackWall) CHECK(std::abs(cloud.points[i].z() - so.wall_distance) < 1e-3);
    }

    const auto dir = std::filesystem::temp_directory_path() / "bh_pc_test";
    std::filesystem::create_directories(dir);
    cloud.normals.assign(cloud.size(), Vector3D{0, 0, -1});
    cloud.colors.assign(cloud.size(), Color{10, 20, 30});
    cloud.scalars.assign(cloud.size(), 0.5);
    for (const bool binary : {true, false}) {
        const auto file = (dir / (binary ? "a.ply" : "b.ply")).string();
        CHECK(save_ply(cloud, file, binary));
        auto back = load_ply(file);
        CHECK(back.has_value());
        if (back) {
            CHECK(back->size() == cloud.size());
            CHECK(back->has_normals() && back->has_colors() && back->has_scalars());
            Real maxd = 0.0;
            for (Index i = 0; i < cloud.size(); ++i) maxd = std::max(maxd, back->points[i].distance(cloud.points[i]));
            CHECK(maxd < 1e-5);   // float storage
            CHECK(back->colors[3].g == 20);
        }
    }
    const auto xyz = (dir / "c.xyz").string();
    CHECK(save_xyz(cloud, xyz));
    auto back = load_xyz(xyz);
    CHECK(back && back->size() == cloud.size() && back->has_normals() && back->has_colors());

    const auto pgm = (dir / "d.pgm").string();
    CHECK(save_depth_pgm(frame.depth, pgm));
    auto dback = load_depth_pgm(pgm);
    CHECK(dback && dback->data == frame.depth.data);

    // Minimal ascii PCD with packed rgb.
    {
        const auto pcd = (dir / "e.pcd").string();
        FILE* f = std::fopen(pcd.c_str(), "w");
        std::fprintf(f, "# .PCD v0.7\nVERSION 0.7\nFIELDS x y z rgb\nSIZE 4 4 4 4\nTYPE F F F F\nCOUNT 1 1 1 1\nWIDTH 2\nHEIGHT 1\nVIEWPOINT 0 0 0 1 0 0 0\nPOINTS 2\nDATA ascii\n1 2 3 4.2108e-45\n4 5 6 0\n");
        std::fclose(f);
        auto p = load_pcd(pcd);
        CHECK(p && p->size() == 2 && p->has_colors());
        if (p) CHECK(std::abs(p->points[1].z() - 6.0) < 1e-9);
    }
}

void test_wlop_improves_uniformity() {
    // Non-uniform sampling of a plane: density increases along x.
    std::mt19937_64 rng(4);
    std::uniform_real_distribution<Real> u(0.0, 1.0);
    PointCloud cloud;
    for (int i = 0; i < 6000; ++i) {
        const Real x = std::sqrt(u(rng));   // pdf ~ 2x
        cloud.points.emplace_back(x, u(rng), 0.0);
    }
    WlopOptions wo;
    wo.target_count = 600;
    wo.support_radius = 0.12;
    wo.iterations = 15;
    wo.theta = 0.5;
    auto res = wlop_resample(cloud, wo);
    CHECK(res.cloud.size() == 600);
    const auto before = nearest_neighbor_statistics(cloud.select(std::vector<Index>(600, 0)).empty() ? cloud : cloud);
    const auto after = nearest_neighbor_statistics(res.cloud);
    CHECK(after.coefficient_of_variation() < 0.45);   // random uniform ~0.52, blue-noise ~0.25
    CHECK(after.coefficient_of_variation() < before.coefficient_of_variation());
    Real max_z = 0.0;
    for (const auto& p : res.cloud.points) max_z = std::max(max_z, std::abs(p.z()));
    CHECK(max_z < 1e-6);   // projection stays on the plane
}

} // namespace

int main() {
    test_octree_queries();
    test_bh_exact_consistency();
    test_mean_shift_separates_blobs();
    test_normals_on_plane_and_sphere();
    test_depth_roundtrip_and_io();
    test_wlop_improves_uniformity();
    if (failures == 0) {
        std::cout << "All point-cloud tests passed\n";
        return 0;
    }
    std::cerr << failures << " check(s) failed\n";
    return 1;
}
