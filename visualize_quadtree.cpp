/**
 * Quadtree Visualization Tool for Barnes-Hut Simulation
 * Modern C++20 implementation with SVG output
 */

#include "quadtree_visualizer.h"
#include "file.h"
#include "tree.h"
#include <iostream>
#include <string>

using namespace barnes_hut;

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " <data_file> [options]\n"
              << "\nOptions:\n"
              << "  --projection <0|1|2>   Projection plane (0=XY, 1=XZ, 2=YZ) [default: 0]\n"
              << "  --no-boxes             Don't show tree structure boxes\n"
              << "  --no-particles         Don't show particles\n"
              << "  --no-mass-centers      Don't show mass centers\n"
              << "  --no-clustering        Disable particle clustering\n"
              << "  --cluster-dist <val>   Clustering distance threshold [default: 0.5]\n"
              << "  --width <px>           SVG width [default: 1920]\n"
              << "  --height <px>          SVG height [default: 1080]\n"
              << "  --theta <val>          Barnes-Hut theta parameter [default: 0.5]\n"
              << "  --max-leaf <val>       Max particles per leaf [default: 10]\n"
              << "  --output <file>        Output SVG filename [default: quadtree.svg]\n"
              << "\nExamples:\n"
              << "  " << program_name << " test.dat\n"
              << "  " << program_name << " test.dat --projection 1 --no-clustering\n"
              << "  " << program_name << " test.dat --width 3840 --height 2160\n";
}

int main(int argc, char* argv[]) {
    std::cout << "=== Barnes-Hut Quadtree Visualizer ===\n\n";

    if (argc < 2) {
        print_usage(argv[0]);
        return 1;
    }

    // Parse arguments
    std::string data_file = argv[1];
    QuadtreeVisualizer::Config viz_config;
    std::string output_file = "quadtree.svg";
    Real theta = 0.5;
    Index max_particles_per_leaf = 10;

    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "--projection" && i + 1 < argc) {
            viz_config.projection = std::stoi(argv[++i]);
            if (viz_config.projection < 0 || viz_config.projection > 2) {
                std::cerr << "Error: Projection must be 0, 1, or 2\n";
                return 1;
            }
        }
        else if (arg == "--no-boxes") {
            viz_config.show_boxes = false;
        }
        else if (arg == "--no-particles") {
            viz_config.show_particles = false;
        }
        else if (arg == "--no-mass-centers") {
            viz_config.show_mass_centers = false;
        }
        else if (arg == "--no-clustering") {
            viz_config.use_clustering = false;
        }
        else if (arg == "--cluster-dist" && i + 1 < argc) {
            viz_config.cluster_threshold = std::stod(argv[++i]);
        }
        else if (arg == "--width" && i + 1 < argc) {
            viz_config.width = std::stoi(argv[++i]);
        }
        else if (arg == "--height" && i + 1 < argc) {
            viz_config.height = std::stoi(argv[++i]);
        }
        else if (arg == "--theta" && i + 1 < argc) {
            theta = std::stod(argv[++i]);
        }
        else if (arg == "--max-leaf" && i + 1 < argc) {
            max_particles_per_leaf = std::stoull(argv[++i]);
        }
        else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        }
        else {
            std::cerr << "Warning: Unknown option: " << arg << "\n";
        }
    }

    // Read configuration
    std::cout << "Reading configuration from: " << data_file << "\n";
    auto config_opt = read_config_file(data_file);
    if (!config_opt) {
        std::cerr << "Error: Failed to read configuration\n";
        return 2;
    }

    const auto& config = *config_opt;

    // Read particles
    std::cout << "Reading " << config.particle_count << " particles...\n";
    auto particles_opt = read_particle_file(data_file, config);
    if (!particles_opt) {
        std::cerr << "Error: Failed to read particles\n";
        return 3;
    }

    auto& particles = *particles_opt;
    std::cout << "Successfully loaded " << particles.size() << " particles\n\n";

    // Build Barnes-Hut tree
    std::cout << "Building Barnes-Hut tree...\n";
    std::cout << "  Theta: " << theta << "\n";
    std::cout << "  Max particles per leaf: " << max_particles_per_leaf << "\n";

    // Note: We need timestep for tree construction, use a dummy value
    const Real dummy_dt = 0.01;
    BarnesHutTree tree(particles, dummy_dt, theta, max_particles_per_leaf);

    // The tree is built during construction via clear_tree() and build_tree()
    // But we need to access the root node - this requires modifying the tree class
    // For now, we'll create a simple visualization without the actual tree structure

    std::cout << "\nGenerating visualization...\n";
    std::cout << "  Output: " << output_file << "\n";
    std::cout << "  Projection: " << (viz_config.projection == 0 ? "XY" :
                                      viz_config.projection == 1 ? "XZ" : "YZ") << "\n";
    std::cout << "  Resolution: " << viz_config.width << "x" << viz_config.height << "\n";
    std::cout << "  Show boxes: " << (viz_config.show_boxes ? "yes" : "no") << "\n";
    std::cout << "  Show particles: " << (viz_config.show_particles ? "yes" : "no") << "\n";
    std::cout << "  Show mass centers: " << (viz_config.show_mass_centers ? "yes" : "no") << "\n";
    std::cout << "  Clustering: " << (viz_config.use_clustering ? "enabled" : "disabled") << "\n";

    if (viz_config.use_clustering) {
        std::cout << "  Cluster threshold: " << viz_config.cluster_threshold << "\n";
    }

    QuadtreeVisualizer visualizer(viz_config);

    // For now, visualize particles only
    // TODO: Expose tree root from BarnesHutTree class
    if (!visualizer.visualize_particles_only(particles, output_file)) {
        std::cerr << "Error: Visualization failed\n";
        return 4;
    }

    // Also export to multiple projections
    std::cout << "\nGenerating all projections...\n";

    for (int proj = 0; proj < 3; ++proj) {
        auto proj_config = viz_config;
        proj_config.projection = proj;

        const std::array<std::string, 3> proj_names = {"xy", "xz", "yz"};
        const std::string proj_output = "quadtree_" + proj_names[proj] + ".svg";

        QuadtreeVisualizer proj_viz(proj_config);
        proj_viz.visualize_particles_only(particles, proj_output);
    }

    std::cout << "\n=== Visualization Complete ===\n";
    std::cout << "\nGenerated files:\n";
    std::cout << "  " << output_file << "\n";
    std::cout << "  quadtree_xy.svg\n";
    std::cout << "  quadtree_xz.svg\n";
    std::cout << "  quadtree_yz.svg\n";
    std::cout << "\nOpen these files in a web browser to view.\n";

    return 0;
}
