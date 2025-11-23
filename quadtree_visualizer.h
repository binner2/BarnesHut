#pragma once

#include "tree.h"
#include "particle.h"
#include <string>
#include <fstream>
#include <array>
#include <cmath>

namespace barnes_hut {

/**
 * Modern C++20 Quadtree/Octree Visualizer
 *
 * Features:
 * - SVG output for scalable visualization
 * - Color-coded tree levels
 * - Clustered particle rendering
 * - Mass center visualization
 * - Transparency for overlapping regions
 */
class QuadtreeVisualizer {
public:
    // Color scheme for different tree levels
    struct ColorScheme {
        std::array<std::string, 10> level_colors = {
            "#FF6B6B",  // Level 0 - Red
            "#4ECDC4",  // Level 1 - Teal
            "#45B7D1",  // Level 2 - Blue
            "#96CEB4",  // Level 3 - Green
            "#FFEAA7",  // Level 4 - Yellow
            "#DFE6E9",  // Level 5 - Light Gray
            "#FF7675",  // Level 6 - Light Red
            "#74B9FF",  // Level 7 - Light Blue
            "#A29BFE",  // Level 8 - Purple
            "#FD79A8"   // Level 9 - Pink
        };

        std::string particle_color = "#FFFFFF";      // White
        std::string cluster_color = "#FFA500";       // Orange
        std::string mass_center_color = "#00FF00";   // Green
        std::string background_color = "#1A1A2E";    // Dark Blue
    };

    // Visualization configuration
    struct Config {
        int width = 1920;                  // SVG width
        int height = 1080;                 // SVG height
        int margin = 50;                   // Margin around plot
        double min_particle_size = 2.0;    // Minimum particle radius
        double max_particle_size = 10.0;   // Maximum particle radius
        double box_opacity = 0.3;          // Tree box transparency
        double particle_opacity = 0.8;     // Particle transparency
        bool show_particles = true;        // Show individual particles
        bool show_boxes = true;            // Show tree boxes
        bool show_mass_centers = true;     // Show mass centers
        bool use_clustering = true;        // Cluster nearby particles
        double cluster_threshold = 0.5;    // Distance for clustering
        int projection = 0;                // 0=XY, 1=XZ, 2=YZ
    };

    QuadtreeVisualizer(const Config& config = Config{})
        : config_(config), colors_(ColorScheme{}) {}

    // Main visualization function
    bool visualize_tree(
        const Node* root,
        std::span<const Particle> particles,
        std::string_view filename) const;

    // Visualize only particles (no tree structure)
    bool visualize_particles_only(
        std::span<const Particle> particles,
        std::string_view filename) const;

    // Visualize with clustering
    bool visualize_clustered(
        const Node* root,
        std::span<const Particle> particles,
        std::string_view filename) const;

    // Export tree data to CSV for analysis
    bool export_tree_data(
        const Node* root,
        std::string_view filename) const;

private:
    Config config_;
    ColorScheme colors_;

    // Helper: Project 3D coordinates to 2D
    struct Point2D {
        double x, y;
    };

    [[nodiscard]] Point2D project(const Vector3D& pos) const noexcept {
        switch (config_.projection) {
            case 1: return {pos[0], pos[2]}; // XZ
            case 2: return {pos[1], pos[2]}; // YZ
            default: return {pos[0], pos[1]}; // XY
        }
    }

    // Helper: Find bounding box for normalization
    struct BoundingBox {
        double min_x, max_x, min_y, max_y;
    };

    [[nodiscard]] BoundingBox find_bounds(std::span<const Particle> particles) const {
        BoundingBox bbox{1e10, -1e10, 1e10, -1e10};

        for (const auto& p : particles) {
            auto [x, y] = project(p.position());
            bbox.min_x = std::min(bbox.min_x, x);
            bbox.max_x = std::max(bbox.max_x, x);
            bbox.min_y = std::min(bbox.min_y, y);
            bbox.max_y = std::max(bbox.max_y, y);
        }

        return bbox;
    }

    // Helper: Transform world coordinates to SVG coordinates
    [[nodiscard]] Point2D world_to_svg(
        double x, double y,
        const BoundingBox& bbox) const noexcept {

        const double plot_width = config_.width - 2 * config_.margin;
        const double plot_height = config_.height - 2 * config_.margin;

        const double scale_x = plot_width / (bbox.max_x - bbox.min_x);
        const double scale_y = plot_height / (bbox.max_y - bbox.min_y);
        const double scale = std::min(scale_x, scale_y);

        const double svg_x = config_.margin + (x - bbox.min_x) * scale;
        const double svg_y = config_.height - config_.margin - (y - bbox.min_y) * scale;

        return {svg_x, svg_y};
    }

    // Helper: Calculate particle size based on mass
    [[nodiscard]] double calculate_particle_size(double mass, double min_mass, double max_mass) const noexcept {
        if (max_mass == min_mass) return config_.min_particle_size;

        const double normalized = (mass - min_mass) / (max_mass - min_mass);
        return config_.min_particle_size +
               normalized * (config_.max_particle_size - config_.min_particle_size);
    }

    // SVG writing helpers
    void write_svg_header(std::ofstream& file) const {
        file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
             << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
             << "width=\"" << config_.width << "\" "
             << "height=\"" << config_.height << "\" "
             << "viewBox=\"0 0 " << config_.width << " " << config_.height << "\">\n";

        // Background
        file << "  <rect width=\"" << config_.width << "\" "
             << "height=\"" << config_.height << "\" "
             << "fill=\"" << colors_.background_color << "\"/>\n";
    }

    void write_svg_footer(std::ofstream& file) const {
        file << "</svg>\n";
    }

    // Recursive tree box drawing
    void draw_node_boxes(
        std::ofstream& file,
        const Node* node,
        const BoundingBox& bbox) const {

        if (!node || node->type == NodeType::Empty) return;

        // Get node bounds in 2D projection
        const auto center_2d = project(node->geo_center);
        const double half_size = node->size / 2.0;

        Vector3D min_corner = node->geo_center;
        Vector3D max_corner = node->geo_center;

        for (int i = 0; i < NDIM; ++i) {
            min_corner[i] -= half_size;
            max_corner[i] += half_size;
        }

        const auto min_2d = project(min_corner);
        const auto max_2d = project(max_corner);

        const auto svg_min = world_to_svg(min_2d.x, min_2d.y, bbox);
        const auto svg_max = world_to_svg(max_2d.x, max_2d.y, bbox);

        const double width = std::abs(svg_max.x - svg_min.x);
        const double height = std::abs(svg_min.y - svg_max.y);

        // Select color based on level
        const auto& color = colors_.level_colors[node->level % colors_.level_colors.size()];

        // Draw rectangle
        file << "  <rect x=\"" << std::min(svg_min.x, svg_max.x) << "\" "
             << "y=\"" << std::min(svg_min.y, svg_max.y) << "\" "
             << "width=\"" << width << "\" "
             << "height=\"" << height << "\" "
             << "fill=\"none\" "
             << "stroke=\"" << color << "\" "
             << "stroke-width=\"1\" "
             << "opacity=\"" << config_.box_opacity << "\"/>\n";

        // Recursively draw children
        for (const auto& child : node->children) {
            if (child) {
                draw_node_boxes(file, child.get(), bbox);
            }
        }
    }

    // Draw particles with clustering
    void draw_particles_clustered(
        std::ofstream& file,
        std::span<const Particle> particles,
        const BoundingBox& bbox) const {

        // Find mass range for size scaling
        double min_mass = 1e10, max_mass = -1e10;
        for (const auto& p : particles) {
            min_mass = std::min(min_mass, p.mass());
            max_mass = std::max(max_mass, p.mass());
        }

        // Track which particles are clustered
        std::vector<bool> clustered(particles.size(), false);

        for (size_t i = 0; i < particles.size(); ++i) {
            if (clustered[i]) continue;

            const auto& p1 = particles[i];
            const auto pos1_2d = project(p1.position());

            // Find nearby particles for clustering
            std::vector<size_t> cluster_members{i};
            double total_mass = p1.mass();
            Vector3D cluster_center = p1.position();

            if (config_.use_clustering) {
                for (size_t j = i + 1; j < particles.size(); ++j) {
                    if (clustered[j]) continue;

                    const auto& p2 = particles[j];
                    const double dist = (p1.position() - p2.position()).length();

                    if (dist < config_.cluster_threshold) {
                        cluster_members.push_back(j);
                        clustered[j] = true;
                        total_mass += p2.mass();
                        cluster_center += p2.position();
                    }
                }

                cluster_center = cluster_center / static_cast<double>(cluster_members.size());
            }

            // Draw cluster or single particle
            const auto center_2d = project(cluster_center);
            const auto svg_pos = world_to_svg(center_2d.x, center_2d.y, bbox);

            const double radius = calculate_particle_size(total_mass, min_mass, max_mass);
            const auto& color = cluster_members.size() > 1 ?
                              colors_.cluster_color : colors_.particle_color;

            file << "  <circle cx=\"" << svg_pos.x << "\" "
                 << "cy=\"" << svg_pos.y << "\" "
                 << "r=\"" << radius << "\" "
                 << "fill=\"" << color << "\" "
                 << "opacity=\"" << config_.particle_opacity << "\"/>\n";

            // Add label for large clusters
            if (cluster_members.size() > 5) {
                file << "  <text x=\"" << svg_pos.x << "\" "
                     << "y=\"" << svg_pos.y + 5 << "\" "
                     << "font-size=\"10\" "
                     << "fill=\"white\" "
                     << "text-anchor=\"middle\">"
                     << cluster_members.size() << "</text>\n";
            }
        }
    }

    // Draw mass centers
    void draw_mass_centers(
        std::ofstream& file,
        const Node* node,
        const BoundingBox& bbox) const {

        if (!node || node->type == NodeType::Empty) return;
        if (node->type == NodeType::Leaf) return; // Skip leaf nodes

        const auto center_2d = project(node->mass_center);
        const auto svg_pos = world_to_svg(center_2d.x, center_2d.y, bbox);

        // Draw mass center as a cross
        const double size = 5.0;
        file << "  <line x1=\"" << svg_pos.x - size << "\" "
             << "y1=\"" << svg_pos.y << "\" "
             << "x2=\"" << svg_pos.x + size << "\" "
             << "y2=\"" << svg_pos.y << "\" "
             << "stroke=\"" << colors_.mass_center_color << "\" "
             << "stroke-width=\"2\"/>\n";
        file << "  <line x1=\"" << svg_pos.x << "\" "
             << "y1=\"" << svg_pos.y - size << "\" "
             << "x2=\"" << svg_pos.x << "\" "
             << "y2=\"" << svg_pos.y + size << "\" "
             << "stroke=\"" << colors_.mass_center_color << "\" "
             << "stroke-width=\"2\"/>\n";

        // Recursively draw children mass centers
        for (const auto& child : node->children) {
            if (child) {
                draw_mass_centers(file, child.get(), bbox);
            }
        }
    }
};

} // namespace barnes_hut
