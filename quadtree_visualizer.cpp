#include "quadtree_visualizer.h"
#include <iostream>
#include <iomanip>
#include <sstream>

namespace barnes_hut {

bool QuadtreeVisualizer::visualize_tree(
    const Node* root,
    std::span<const Particle> particles,
    std::string_view filename) const {

    std::ofstream file(filename.data());
    if (!file) {
        std::cerr << "Error: Could not create visualization file: " << filename << "\n";
        return false;
    }

    const auto bbox = find_bounds(particles);

    write_svg_header(file);

    // Add title
    const std::array<std::string, 3> proj_names = {"XY", "XZ", "YZ"};
    file << "  <text x=\"" << config_.width / 2 << "\" y=\"30\" "
         << "font-size=\"24\" fill=\"white\" text-anchor=\"middle\">"
         << "Barnes-Hut Quadtree (" << proj_names[config_.projection] << " Projection)"
         << "</text>\n";

    // Draw tree boxes
    if (config_.show_boxes && root) {
        file << "  <!-- Tree Structure -->\n";
        draw_node_boxes(file, root, bbox);
    }

    // Draw mass centers
    if (config_.show_mass_centers && root) {
        file << "  <!-- Mass Centers -->\n";
        draw_mass_centers(file, root, bbox);
    }

    // Draw particles
    if (config_.show_particles) {
        file << "  <!-- Particles -->\n";
        draw_particles_clustered(file, particles, bbox);
    }

    // Add legend
    file << "  <!-- Legend -->\n";
    const int legend_x = config_.width - 200;
    int legend_y = 80;

    file << "  <rect x=\"" << legend_x - 10 << "\" y=\"" << legend_y - 25 << "\" "
         << "width=\"180\" height=\"150\" "
         << "fill=\"#000000\" opacity=\"0.7\" rx=\"5\"/>\n";

    file << "  <text x=\"" << legend_x << "\" y=\"" << legend_y << "\" "
         << "font-size=\"14\" fill=\"white\" font-weight=\"bold\">Legend</text>\n";
    legend_y += 25;

    // Particle legend
    file << "  <circle cx=\"" << legend_x << "\" cy=\"" << legend_y << "\" "
         << "r=\"5\" fill=\"" << colors_.particle_color << "\"/>\n";
    file << "  <text x=\"" << legend_x + 15 << "\" y=\"" << legend_y + 5 << "\" "
         << "font-size=\"12\" fill=\"white\">Particles</text>\n";
    legend_y += 20;

    // Cluster legend
    if (config_.use_clustering) {
        file << "  <circle cx=\"" << legend_x << "\" cy=\"" << legend_y << "\" "
             << "r=\"5\" fill=\"" << colors_.cluster_color << "\"/>\n";
        file << "  <text x=\"" << legend_x + 15 << "\" y=\"" << legend_y + 5 << "\" "
             << "font-size=\"12\" fill=\"white\">Clusters</text>\n";
        legend_y += 20;
    }

    // Mass center legend
    if (config_.show_mass_centers) {
        file << "  <line x1=\"" << legend_x - 5 << "\" y1=\"" << legend_y << "\" "
             << "x2=\"" << legend_x + 5 << "\" y2=\"" << legend_y << "\" "
             << "stroke=\"" << colors_.mass_center_color << "\" stroke-width=\"2\"/>\n";
        file << "  <line x1=\"" << legend_x << "\" y1=\"" << legend_y - 5 << "\" "
             << "x2=\"" << legend_x << "\" y2=\"" << legend_y + 5 << "\" "
             << "stroke=\"" << colors_.mass_center_color << "\" stroke-width=\"2\"/>\n";
        file << "  <text x=\"" << legend_x + 15 << "\" y=\"" << legend_y + 5 << "\" "
             << "font-size=\"12\" fill=\"white\">Mass Centers</text>\n";
        legend_y += 20;
    }

    // Tree levels legend
    if (config_.show_boxes) {
        file << "  <text x=\"" << legend_x << "\" y=\"" << legend_y << "\" "
             << "font-size=\"12\" fill=\"white\">Tree Levels:</text>\n";
        legend_y += 15;

        for (int i = 0; i < 3 && i < static_cast<int>(colors_.level_colors.size()); ++i) {
            file << "  <rect x=\"" << legend_x << "\" y=\"" << legend_y - 8 << "\" "
                 << "width=\"15\" height=\"10\" "
                 << "fill=\"none\" stroke=\"" << colors_.level_colors[i] << "\" "
                 << "stroke-width=\"2\"/>\n";
            file << "  <text x=\"" << legend_x + 20 << "\" y=\"" << legend_y << "\" "
                 << "font-size=\"10\" fill=\"white\">Level " << i << "</text>\n";
            legend_y += 15;
        }
    }

    // Statistics
    file << "  <!-- Statistics -->\n";
    file << "  <text x=\"20\" y=\"" << config_.height - 20 << "\" "
         << "font-size=\"12\" fill=\"#AAAAAA\">"
         << "Particles: " << particles.size()
         << "</text>\n";

    write_svg_footer(file);

    std::cout << "Quadtree visualization saved to: " << filename << "\n";
    return true;
}

bool QuadtreeVisualizer::visualize_particles_only(
    std::span<const Particle> particles,
    std::string_view filename) const {

    auto modified_config = config_;
    modified_config.show_boxes = false;
    modified_config.show_mass_centers = false;

    QuadtreeVisualizer viz(modified_config);
    return viz.visualize_tree(nullptr, particles, filename);
}

bool QuadtreeVisualizer::visualize_clustered(
    const Node* root,
    std::span<const Particle> particles,
    std::string_view filename) const {

    auto modified_config = config_;
    modified_config.use_clustering = true;
    modified_config.cluster_threshold = 0.5;

    QuadtreeVisualizer viz(modified_config);
    return viz.visualize_tree(root, particles, filename);
}

bool QuadtreeVisualizer::export_tree_data(
    const Node* root,
    std::string_view filename) const {

    std::ofstream file(filename.data());
    if (!file) {
        std::cerr << "Error: Could not create export file: " << filename << "\n";
        return false;
    }

    // CSV header
    file << "Level,Type,CenterX,CenterY,CenterZ,Size,Mass,MassCenterX,MassCenterY,MassCenterZ,ParticleCount\n";

    // Recursive lambda to export node data
    std::function<void(const Node*)> export_node = [&](const Node* node) {
        if (!node || node->type == NodeType::Empty) return;

        // Export current node
        file << node->level << ","
             << static_cast<int>(node->type) << ","
             << std::fixed << std::setprecision(6)
             << node->geo_center[0] << ","
             << node->geo_center[1] << ","
             << node->geo_center[2] << ","
             << node->size << ","
             << node->mass << ","
             << node->mass_center[0] << ","
             << node->mass_center[1] << ","
             << node->mass_center[2] << ","
             << node->particle_count << "\n";

        // Export children
        for (const auto& child : node->children) {
            if (child) {
                export_node(child.get());
            }
        }
    };

    export_node(root);

    std::cout << "Tree data exported to: " << filename << "\n";
    return true;
}

} // namespace barnes_hut
