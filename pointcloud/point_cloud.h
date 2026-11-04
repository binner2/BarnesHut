#pragma once
/**
 * Point cloud container, depth-camera back-projection and file I/O.
 *
 * Reuses the repository's Vector3D / Real / Index types so that the
 * Barnes-Hut machinery and the point-cloud pipeline share one numeric core.
 */
#include "../vektor.h"
#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace barnes_hut::pc {

struct Color {
    std::uint8_t r = 0, g = 0, b = 0;
};

struct PointCloud {
    std::vector<Vector3D> points;
    std::vector<Vector3D> normals;   // empty, or one per point
    std::vector<Color> colors;       // empty, or one per point
    std::vector<Real> scalars;       // empty, or one per point (density, label, ...)

    [[nodiscard]] Index size() const noexcept { return points.size(); }
    [[nodiscard]] bool empty() const noexcept { return points.empty(); }
    [[nodiscard]] bool has_normals() const noexcept { return !points.empty() && normals.size() == points.size(); }
    [[nodiscard]] bool has_colors() const noexcept { return !points.empty() && colors.size() == points.size(); }
    [[nodiscard]] bool has_scalars() const noexcept { return !points.empty() && scalars.size() == points.size(); }

    void reserve(Index n);
    void clear();

    // Keep only the listed indices (attributes follow).
    [[nodiscard]] PointCloud select(std::span<const Index> indices) const;

    // Axis-aligned bounding box. Returns false for an empty cloud.
    bool bounding_box(Vector3D& lo, Vector3D& hi) const noexcept;
    [[nodiscard]] Vector3D centroid() const noexcept;
};

// Pinhole model. depth_scale converts raw depth units to metres
// (16-bit millimetre images: 1000; TUM RGB-D PNGs: 5000).
struct CameraIntrinsics {
    Real fx = 525.0;
    Real fy = 525.0;
    Real cx = 319.5;
    Real cy = 239.5;
    int width = 640;
    int height = 480;
    Real depth_scale = 1000.0;

    [[nodiscard]] Vector3D back_project(Real u, Real v, Real z_metres) const noexcept {
        return Vector3D{(u - cx) * z_metres / fx, (v - cy) * z_metres / fy, z_metres};
    }
    // Returns false if the point is behind the camera.
    bool project(const Vector3D& p, Real& u, Real& v) const noexcept {
        if (p.z() <= 0.0) return false;
        u = fx * p.x() / p.z() + cx;
        v = fy * p.y() / p.z() + cy;
        return true;
    }
};

struct DepthImage {
    int width = 0;
    int height = 0;
    std::vector<std::uint16_t> data;   // row-major, raw units

    [[nodiscard]] std::uint16_t at(int u, int v) const noexcept {
        return data[static_cast<std::size_t>(v) * width + u];
    }
    [[nodiscard]] std::uint16_t& at(int u, int v) noexcept {
        return data[static_cast<std::size_t>(v) * width + u];
    }
};

struct DepthToCloudOptions {
    Real min_depth = 0.05;   // metres
    Real max_depth = 10.0;   // metres
    int stride = 1;          // pixel subsampling step
};

// Back-project a depth image. If pixel_of_point is given it receives (u, v)
// for every produced point, which lets callers map results back to the image.
PointCloud depth_to_point_cloud(const DepthImage& depth,
                                const CameraIntrinsics& K,
                                const DepthToCloudOptions& opt = {},
                                std::vector<std::pair<int, int>>* pixel_of_point = nullptr);

// --- Depth image readers ---------------------------------------------------
std::optional<DepthImage> load_depth_pgm(std::string_view filename);           // P2 / P5, 8 or 16 bit
std::optional<DepthImage> load_depth_raw(std::string_view filename, int width, int height); // uint16 LE
std::optional<DepthImage> load_depth_png(std::string_view filename);           // needs BH_HAVE_PNG
std::optional<DepthImage> load_depth_image(std::string_view filename);         // by extension
bool save_depth_pgm(const DepthImage& depth, std::string_view filename);
[[nodiscard]] bool png_support_compiled() noexcept;

// --- Intrinsics -------------------------------------------------------------
// Text file: "fx fy cx cy [width height [depth_scale]]", '#' comments allowed.
std::optional<CameraIntrinsics> load_intrinsics(std::string_view filename);

// --- Point cloud readers / writers ------------------------------------------
std::optional<PointCloud> load_point_cloud(std::string_view filename);        // by extension
std::optional<PointCloud> load_xyz(std::string_view filename);                // x y z [nx ny nz] [r g b]
std::optional<PointCloud> load_ply(std::string_view filename);                // ascii / binary_little_endian
std::optional<PointCloud> load_pcd(std::string_view filename);                // ascii / binary
bool save_point_cloud(const PointCloud& cloud, std::string_view filename);    // by extension
bool save_xyz(const PointCloud& cloud, std::string_view filename);
bool save_ply(const PointCloud& cloud, std::string_view filename, bool binary = true);

// Map a scalar field to a colour ramp (viridis-like) for visual inspection.
void colorize_by_scalar(PointCloud& cloud);
// Distinct colours for integer labels (negative label = noise, grey).
void colorize_by_label(PointCloud& cloud);

} // namespace barnes_hut::pc
