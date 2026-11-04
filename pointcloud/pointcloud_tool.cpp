/**
 * Command-line front end for the depth-camera point-cloud pipeline.
 *
 *   pointcloud_tool <command> [options]
 *
 * Commands: info, synthetic, depth2cloud, downsample, normals, sor, density,
 *           meanshift, resample, bench. Run without arguments for help.
 */
#include "algorithms.h"
#include "point_cloud.h"
#include "point_octree.h"
#include "synthetic.h"

#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace barnes_hut;
using namespace barnes_hut::pc;

namespace {

struct Args {
    std::string command;
    std::vector<std::string> positional;
    std::map<std::string, std::string> options;

    [[nodiscard]] bool has(const std::string& k) const { return options.count(k) > 0; }
    [[nodiscard]] std::string get(const std::string& k, const std::string& def = "") const {
        auto it = options.find(k);
        return it == options.end() ? def : it->second;
    }
    [[nodiscard]] double num(const std::string& k, double def) const {
        auto it = options.find(k);
        return it == options.end() ? def : std::stod(it->second);
    }
    [[nodiscard]] long integer(const std::string& k, long def) const {
        auto it = options.find(k);
        return it == options.end() ? def : std::stol(it->second);
    }
};

Args parse(int argc, char** argv) {
    Args a;
    if (argc >= 2) a.command = argv[1];
    for (int i = 2; i < argc; ++i) {
        std::string s = argv[i];
        if (s.rfind("--", 0) == 0) {
            const std::string key = s.substr(2);
            auto is_value = [](const std::string& t) {
                if (t.empty()) return false;
                if (t[0] != '-') return true;
                return t.size() > 1 && (std::isdigit(static_cast<unsigned char>(t[1])) || t[1] == '.');   // negative number
            };
            if (i + 1 < argc && is_value(argv[i + 1])) a.options[key] = argv[++i];
            else a.options[key] = "1";
        } else if (s == "-o" && i + 1 < argc) {
            a.options["output"] = argv[++i];
        } else {
            a.positional.push_back(s);
        }
    }
    return a;
}

void usage() {
    std::cout <<
        "pointcloud_tool - Barnes-Hut accelerated depth-camera point-cloud pipeline\n\n"
        "  info        <cloud>                               statistics, octree summary\n"
        "  synthetic   -o out.ply [--pgm depth.pgm] [--noise c] [--labels]\n"
        "                                                     render a synthetic depth frame (640x480)\n"
        "  depth2cloud <depth.png|.pgm|.raw> -o out.ply [--intrinsics file | --fx --fy --cx --cy]\n"
        "              [--scale 1000] [--width W --height H] [--stride 1] [--min 0.05] [--max 10]\n"
        "  downsample  <cloud> --voxel 0.01 -o out.ply\n"
        "  normals     <cloud> [--k 20 | --radius r] [--viewpoint x,y,z] -o out.ply\n"
        "  sor         <cloud> [--k 20] [--std 2.0] -o out.ply        statistical outlier removal\n"
        "  density     <cloud> --h 0.02 [--theta 0.5] [--kappa 0] [--cutoff 4] [--compare] -o out.ply\n"
        "  meanshift   <cloud> --h 0.05 [--theta 0.5] [--kappa 0] [--min-size 50] -o out.ply\n"
        "  resample    <cloud> --count 5000 --h 0.05 [--mu 0.45] [--iters 20] [--theta 0.5] [--kappa 0] -o out.ply\n"
        "  bench       [--sizes 20000,50000,100000] [--thetas 0,0.3,0.5,0.7,1.0] [--kappas 0,1,0.5] [--h 0.03] --csv out.csv\n\n"
        "Common: --leaf N (octree leaf capacity, default 16), --threads T\n"
        "theta: Barnes-Hut opening angle (0 = exact). kappa: kernel-aware limit, approximate only cells\n"
        "narrower than kappa*h (0 = geometric criterion only).\n";
}

std::optional<Vector3D> parse_vec(const std::string& s) {
    if (s.empty()) return std::nullopt;
    std::string t = s;
    for (auto& c : t) if (c == ',') c = ' ';
    std::istringstream is(t);
    Real x, y, z;
    if (is >> x >> y >> z) return Vector3D{x, y, z};
    return std::nullopt;
}

PointOctree::Params tree_params(const Args& a) {
    PointOctree::Params p;
    p.leaf_capacity = static_cast<Index>(a.integer("leaf", 16));
    return p;
}

std::optional<PointCloud> load_or_fail(const Args& a) {
    if (a.positional.empty()) {
        std::cerr << "Error: input cloud missing\n";
        return std::nullopt;
    }
    return load_point_cloud(a.positional.front());
}

bool save_or_fail(const PointCloud& c, const Args& a) {
    const std::string out = a.get("output");
    if (out.empty()) {
        std::cerr << "Error: -o <output> missing\n";
        return false;
    }
    if (!save_point_cloud(c, out)) return false;
    std::cout << "Wrote " << c.size() << " points to " << out << "\n";
    return true;
}

void print_counters(const char* what, const BhCounters& c, Index n) {
    std::cout << what << ": " << std::fixed << std::setprecision(4) << c.seconds << " s"
              << " | far(node) evals " << c.far_evaluations << " | near(point) evals " << c.near_evaluations
              << " | per query " << std::setprecision(1)
              << (n ? static_cast<double>(c.far_evaluations + c.near_evaluations) / static_cast<double>(n) : 0.0) << "\n";
}

// ---------------------------------------------------------------------------

int cmd_info(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    Vector3D lo, hi;
    cloud->bounding_box(lo, hi);
    std::cout << "Points: " << cloud->size() << "  normals: " << (cloud->has_normals() ? "yes" : "no")
              << "  colors: " << (cloud->has_colors() ? "yes" : "no") << "  scalars: " << (cloud->has_scalars() ? "yes" : "no") << "\n"
              << "Bounding box: [" << lo << "] .. [" << hi << "]\n";
    PointOctree tree(cloud->points, tree_params(a));
    const auto& s = tree.statistics();
    std::cout << "Octree: nodes " << s.node_count << ", leaves " << s.leaf_count << ", depth " << s.depth
              << ", root width " << s.root_width << " m, build " << std::setprecision(4) << s.build_seconds << " s\n";
    const auto nn = nearest_neighbor_statistics(*cloud);
    std::cout << "Nearest-neighbour distance: mean " << nn.mean << " m, std " << nn.stddev << " m, CoV "
              << nn.coefficient_of_variation() << "\n";
    return 0;
}

int cmd_synthetic(const Args& a) {
    CameraIntrinsics K;   // 640x480, fx=fy=525, TUM-style defaults, depth_scale 1000 (mm)
    K.depth_scale = a.num("scale", 1000.0);
    SyntheticSceneOptions so;
    so.noise_coefficient = a.num("noise", so.noise_coefficient);
    so.dropout_fraction = a.num("dropout", so.dropout_fraction);
    so.seed = static_cast<unsigned>(a.integer("seed", so.seed));
    const auto frame = render_synthetic_depth(K, so);
    if (a.has("pgm")) {
        save_depth_pgm(frame.depth, a.get("pgm"));
        std::cout << "Wrote depth image " << a.get("pgm") << " (16-bit PGM, raw unit = 1/" << K.depth_scale << " m)\n";
    }
    std::vector<std::pair<int, int>> pix;
    DepthToCloudOptions dopt;
    dopt.stride = static_cast<int>(a.integer("stride", 1));
    PointCloud cloud = depth_to_point_cloud(frame.depth, K, dopt, &pix);
    if (a.has("labels")) {
        cloud.scalars.resize(cloud.size());
        for (Index i = 0; i < cloud.size(); ++i)
            cloud.scalars[i] = frame.surface[static_cast<std::size_t>(pix[i].second) * K.width + pix[i].first];
        colorize_by_label(cloud);
    }
    return save_or_fail(cloud, a) ? 0 : 3;
}

int cmd_depth2cloud(const Args& a) {
    if (a.positional.empty()) {
        std::cerr << "Error: depth image missing\n";
        return 2;
    }
    CameraIntrinsics K;
    if (a.has("intrinsics")) {
        auto k = load_intrinsics(a.get("intrinsics"));
        if (!k) return 2;
        K = *k;
    }
    if (a.has("fx")) K.fx = a.num("fx", K.fx);
    if (a.has("fy")) K.fy = a.num("fy", K.fy);
    if (a.has("cx")) K.cx = a.num("cx", K.cx);
    if (a.has("cy")) K.cy = a.num("cy", K.cy);
    if (a.has("scale")) K.depth_scale = a.num("scale", K.depth_scale);
    const std::string file = a.positional.front();
    std::optional<DepthImage> depth;
    if (file.size() > 4 && file.substr(file.size() - 4) == ".raw")
        depth = load_depth_raw(file, static_cast<int>(a.integer("width", K.width)), static_cast<int>(a.integer("height", K.height)));
    else
        depth = load_depth_image(file);
    if (!depth) return 2;
    K.width = depth->width;
    K.height = depth->height;
    DepthToCloudOptions dopt;
    dopt.stride = static_cast<int>(a.integer("stride", 1));
    dopt.min_depth = a.num("min", dopt.min_depth);
    dopt.max_depth = a.num("max", dopt.max_depth);
    PointCloud cloud = depth_to_point_cloud(*depth, K, dopt);
    std::cout << "Back-projected " << cloud.size() << " of " << depth->width * depth->height << " pixels (fx=" << K.fx
              << ", fy=" << K.fy << ", cx=" << K.cx << ", cy=" << K.cy << ", scale=" << K.depth_scale << ")\n";
    return save_or_fail(cloud, a) ? 0 : 3;
}

int cmd_downsample(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    const Real voxel = a.num("voxel", 0.01);
    Timer t;
    PointCloud out = voxel_downsample(*cloud, voxel);
    std::cout << "Voxel " << voxel << " m: " << cloud->size() << " -> " << out.size() << " points in " << t.elapsed() << " s\n";
    return save_or_fail(out, a) ? 0 : 3;
}

int cmd_normals(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    PointOctree tree(cloud->points, tree_params(a));
    NormalOptions no;
    no.k = static_cast<Index>(a.integer("k", 20));
    no.radius = a.num("radius", 0.0);
    no.viewpoint = a.has("viewpoint") ? parse_vec(a.get("viewpoint")) : std::optional<Vector3D>{Vector3D{0.0}};
    Timer t;
    auto variation = estimate_normals(*cloud, tree, no);
    std::cout << "Normals for " << cloud->size() << " points in " << t.elapsed() << " s (octree build "
              << tree.statistics().build_seconds << " s)\n";
    if (a.has("curvature")) {
        cloud->scalars = std::move(variation);
        colorize_by_scalar(*cloud);
    }
    return save_or_fail(*cloud, a) ? 0 : 3;
}

int cmd_sor(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    PointOctree tree(cloud->points, tree_params(a));
    Timer t;
    auto keep = statistical_outlier_filter(*cloud, tree, static_cast<Index>(a.integer("k", 20)), a.num("std", 2.0));
    std::cout << "SOR kept " << keep.size() << " of " << cloud->size() << " points in " << t.elapsed() << " s\n";
    return save_or_fail(cloud->select(keep), a) ? 0 : 3;
}

int cmd_density(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    PointOctree tree(cloud->points, tree_params(a));
    KernelDensityOptions ko;
    ko.bandwidth = a.num("h", 0.02);
    ko.theta = a.num("theta", 0.5);
    ko.kappa = a.num("kappa", 0.0);
    ko.cutoff_sigmas = a.num("cutoff", 4.0);
    auto res = kernel_density(*cloud, tree, ko);
    print_counters("BH density", res.counters, cloud->size());
    if (a.has("compare")) {
        KernelDensityOptions ex = ko;
        ex.theta = 0.0;
        auto exact = kernel_density(*cloud, tree, ex);
        print_counters("Exact density", exact.counters, cloud->size());
        const auto err = relative_error(res.density, exact.density);
        std::cout << "Relative error vs exact: max " << std::scientific << err.max_relative << ", mean " << err.mean_relative
                  << ", rms " << err.rms_relative << std::defaultfloat << "\n";
    }
    cloud->scalars = std::move(res.density);
    colorize_by_scalar(*cloud);
    return save_or_fail(*cloud, a) ? 0 : 3;
}

int cmd_meanshift(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    PointOctree tree(cloud->points, tree_params(a));
    MeanShiftOptions mo;
    mo.bandwidth = a.num("h", 0.05);
    mo.theta = a.num("theta", 0.5);
    mo.kappa = a.num("kappa", 0.0);
    mo.cutoff_sigmas = a.num("cutoff", 4.0);
    mo.max_iterations = static_cast<int>(a.integer("iters", 50));
    mo.min_cluster_size = static_cast<Index>(a.integer("min-size", 50));
    auto res = mean_shift(*cloud, tree, mo);
    print_counters("BH mean shift", res.counters, cloud->size());
    std::cout << "Clusters: " << res.modes.size() << ", mean iterations per point "
              << static_cast<double>(res.total_iterations) / static_cast<double>(cloud->size()) << "\n";
    for (Index c = 0; c < std::min<Index>(res.modes.size(), 10); ++c)
        std::cout << "  label " << c << ": " << res.cluster_sizes[c] << " points, mode [" << res.modes[c] << "]\n";
    cloud->scalars.resize(cloud->size());
    for (Index i = 0; i < cloud->size(); ++i) cloud->scalars[i] = res.labels[i];
    colorize_by_label(*cloud);
    return save_or_fail(*cloud, a) ? 0 : 3;
}

int cmd_resample(const Args& a) {
    auto cloud = load_or_fail(a);
    if (!cloud) return 2;
    WlopOptions wo;
    wo.target_count = static_cast<Index>(a.integer("count", 0));
    wo.support_radius = a.num("h", 0.05);
    wo.repulsion_mu = a.num("mu", 0.45);
    wo.iterations = static_cast<int>(a.integer("iters", 20));
    wo.theta = a.num("theta", 0.5);
    wo.kappa = a.num("kappa", 0.0);
    wo.seed = static_cast<unsigned>(a.integer("seed", 42));
    const auto before = nearest_neighbor_statistics(*cloud);
    auto res = wlop_resample(*cloud, wo);
    print_counters("BH WLOP", res.counters, res.cloud.size() * static_cast<Index>(wo.iterations));
    const auto after = nearest_neighbor_statistics(res.cloud);
    std::cout << std::fixed << std::setprecision(3) << "NN-distance CoV: input " << before.coefficient_of_variation()
              << " -> resampled " << after.coefficient_of_variation() << " (" << res.cloud.size() << " points)\n";
    return save_or_fail(res.cloud, a) ? 0 : 3;
}

// Benchmark: synthetic frames subsampled to several sizes, theta sweep, exact reference.
int cmd_bench(const Args& a) {
    auto split = [](const std::string& s) {
        std::vector<double> v;
        std::string t = s;
        for (auto& c : t) if (c == ',') c = ' ';
        std::istringstream is(t);
        double x;
        while (is >> x) v.push_back(x);
        return v;
    };
    const auto sizes = split(a.get("sizes", "20000,50000,100000"));
    const auto thetas = split(a.get("thetas", "0,0.3,0.5,0.7,1.0"));
    const auto kappas = split(a.get("kappas", "0,1,0.5"));
    const Real h = a.num("h", 0.03);
    const Real cutoff = a.num("cutoff", 4.0);
    const int repeats = static_cast<int>(a.integer("repeats", 1));
    std::ofstream csv;
    if (a.has("csv")) {
        csv.open(a.get("csv"));
        csv << "algorithm,n,theta,kappa,cutoff_sigmas,leaf,threads,build_s,time_s,far_evals,near_evals,max_rel_err,mean_rel_err,rms_rel_err,extra\n";
    }
    CameraIntrinsics K;
    SyntheticSceneOptions so;
    const auto frame = render_synthetic_depth(K, so);
    PointCloud full = depth_to_point_cloud(frame.depth, K);
    std::cout << "Synthetic frame: " << full.size() << " valid points\n";
    const int threads = [] {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }();

    std::mt19937_64 rng(1);
    for (const double size_d : sizes) {
        const Index n = std::min<Index>(static_cast<Index>(size_d), full.size());
        std::vector<Index> idx(full.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::shuffle(idx.begin(), idx.end(), rng);
        idx.resize(n);
        PointCloud cloud = full.select(idx);
        PointOctree tree(cloud.points, tree_params(a));
        const double build_s = tree.statistics().build_seconds;

        // Reference: exact Gaussian KDE (theta = 0, no cutoff).
        KernelDensityOptions ex;
        ex.bandwidth = h;
        ex.theta = 0.0;
        ex.cutoff_sigmas = 0.0;
        auto exact = kernel_density(cloud, tree, ex);
        std::cout << "n=" << n << " exact KDE " << exact.counters.seconds << " s\n";
        if (csv)
            csv << "kde_exact," << n << ",0,0,0," << a.integer("leaf", 16) << "," << threads << "," << build_s << "," << exact.counters.seconds
                << "," << exact.counters.far_evaluations << "," << exact.counters.near_evaluations << ",0,0,0,\n";

        for (const double theta : thetas) {
            for (const double kappa : kappas) {
                if (theta == 0.0 && kappa != 0.0) continue;   // kappa is irrelevant for the exact sum
                for (const double co : {0.0, cutoff}) {
                    if (theta == 0.0 && co == 0.0) continue;   // identical to the reference
                    if (kappa != 0.0 && co == 0.0) continue;   // keep the sweep small: kappa only with cutoff
                    KernelDensityOptions ko;
                    ko.bandwidth = h;
                    ko.theta = theta;
                    ko.kappa = kappa;
                    ko.cutoff_sigmas = co;
                    double best = 1e300;
                    KernelDensityResult res;
                    for (int r = 0; r < repeats; ++r) {
                        res = kernel_density(cloud, tree, ko);
                        best = std::min(best, res.counters.seconds);
                    }
                    const auto err = relative_error(res.density, exact.density);
                    std::cout << "  KDE theta=" << theta << " kappa=" << kappa << " cutoff=" << co << "h: " << std::fixed
                              << std::setprecision(4) << best << " s, max rel err " << std::scientific << std::setprecision(2)
                              << err.max_relative << ", rms " << err.rms_relative << std::defaultfloat << "\n";
                    if (csv)
                        csv << "kde," << n << "," << theta << "," << kappa << "," << co << "," << a.integer("leaf", 16) << "," << threads
                            << "," << build_s << "," << best << "," << res.counters.far_evaluations << "," << res.counters.near_evaluations
                            << "," << err.max_relative << "," << err.mean_relative << "," << err.rms_relative << ",\n";
                }
            }
        }

        // Mean shift: exact (theta 0, cutoff) vs BH.
        MeanShiftOptions mo;
        mo.bandwidth = a.num("ms-h", 0.08);
        mo.cutoff_sigmas = cutoff;
        mo.min_cluster_size = 50;
        mo.max_iterations = 30;
        mo.theta = 0.0;
        auto ms_exact = mean_shift(cloud, tree, mo);
        if (csv)
            csv << "meanshift," << n << ",0,0," << cutoff << "," << a.integer("leaf", 16) << "," << threads << "," << build_s << ","
                << ms_exact.counters.seconds << "," << ms_exact.counters.far_evaluations << "," << ms_exact.counters.near_evaluations
                << ",0,0,0,clusters=" << ms_exact.modes.size() << "\n";
        std::cout << "  mean shift exact: " << ms_exact.counters.seconds << " s, " << ms_exact.modes.size() << " clusters\n";
        for (const double theta : thetas) {
            if (theta == 0.0) continue;
            mo.theta = theta;
            auto ms = mean_shift(cloud, tree, mo);
            // Label agreement with the exact run (adjusted for permutation via majority mapping).
            std::map<std::pair<int, int>, Index> confusion;
            for (Index i = 0; i < n; ++i) ++confusion[{ms_exact.labels[i], ms.labels[i]}];
            std::map<int, std::pair<int, Index>> best_match;
            for (const auto& [key, cnt] : confusion)
                if (cnt > best_match[key.first].second) best_match[key.first] = {key.second, cnt};
            Index agree = 0;
            for (const auto& [l, m] : best_match) agree += m.second;
            const double agreement = static_cast<double>(agree) / static_cast<double>(n);
            std::cout << "  mean shift theta=" << theta << ": " << ms.counters.seconds << " s, " << ms.modes.size()
                      << " clusters, label agreement " << agreement << "\n";
            if (csv)
                csv << "meanshift," << n << "," << theta << ",0," << cutoff << "," << a.integer("leaf", 16) << "," << threads << ","
                    << build_s << "," << ms.counters.seconds << "," << ms.counters.far_evaluations << "," << ms.counters.near_evaluations
                    << ",0,0,0,clusters=" << ms.modes.size() << ";agreement=" << agreement << "\n";
        }

        // WLOP: exact vs BH, uniformity of the result.
        WlopOptions wo;
        wo.target_count = std::max<Index>(1000, n / 10);
        wo.support_radius = a.num("wlop-h", 0.08);
        wo.iterations = 10;
        for (const double theta : thetas) {
            wo.theta = theta;
            auto w = wlop_resample(cloud, wo);
            const auto nn = nearest_neighbor_statistics(w.cloud);
            std::cout << std::fixed << std::setprecision(4) << "  WLOP theta=" << theta << ": " << w.counters.seconds
                      << " s, NN CoV " << std::setprecision(3) << nn.coefficient_of_variation() << std::defaultfloat << "\n";
            if (csv)
                csv << "wlop," << n << "," << theta << ",0,1," << a.integer("leaf", 16) << "," << threads << "," << build_s << ","
                    << w.counters.seconds << "," << w.counters.far_evaluations << "," << w.counters.near_evaluations
                    << ",0,0,0,nn_cov=" << nn.coefficient_of_variation() << ";m=" << w.cloud.size() << "\n";
        }
    }
    if (csv) std::cout << "Wrote " << a.get("csv") << "\n";
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    const Args a = parse(argc, argv);
#ifdef _OPENMP
    if (a.has("threads")) omp_set_num_threads(static_cast<int>(a.integer("threads", 1)));
#endif
    if (a.command.empty() || a.command == "help" || a.command == "--help") {
        usage();
        return a.command.empty() ? 1 : 0;
    }
    if (a.command == "info") return cmd_info(a);
    if (a.command == "synthetic") return cmd_synthetic(a);
    if (a.command == "depth2cloud") return cmd_depth2cloud(a);
    if (a.command == "downsample") return cmd_downsample(a);
    if (a.command == "normals") return cmd_normals(a);
    if (a.command == "sor") return cmd_sor(a);
    if (a.command == "density") return cmd_density(a);
    if (a.command == "meanshift") return cmd_meanshift(a);
    if (a.command == "resample") return cmd_resample(a);
    if (a.command == "bench") return cmd_bench(a);
    std::cerr << "Unknown command: " << a.command << "\n";
    usage();
    return 1;
}
