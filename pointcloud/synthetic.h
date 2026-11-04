#pragma once
/**
 * Synthetic depth-camera scene for tests and benchmarks: a ray-cast indoor
 * corner (floor, back wall, side wall) with a sphere and a box, observed by a
 * pinhole camera at the origin looking down +z. Depth noise follows the
 * quadratic-in-distance model of structured-light sensors,
 *     sigma_z(z) = noise_coefficient * z^2,
 * followed by quantisation to the 16-bit raw unit of the intrinsics.
 */
#include "point_cloud.h"
#include <cstdint>

namespace barnes_hut::pc {

struct SyntheticSceneOptions {
    Real wall_distance = 3.0;       // back wall z (m)
    Real floor_height = 1.2;        // camera height above floor (m)
    Real side_wall_x = 1.5;         // right side wall at x = +side_wall_x
    Vector3D sphere_center{-0.4, 0.5, 2.0};
    Real sphere_radius = 0.3;
    Vector3D box_min{0.3, 0.2, 1.6};
    Vector3D box_max{0.9, 1.2, 2.2};
    Real noise_coefficient = 0.0019; // sigma_z = c * z^2 (m); ~1.9 mm at 1 m
    Real dropout_fraction = 0.02;    // random invalid pixels (depth = 0)
    unsigned seed = 7;
};

// Ground-truth labels of the scene surfaces, in the order used by the renderer.
enum class SurfaceId : std::uint8_t { None = 0, Floor = 1, BackWall = 2, SideWall = 3, Sphere = 4, Box = 5 };

struct SyntheticDepth {
    DepthImage depth;                       // noisy, quantised
    DepthImage depth_clean;                 // noise-free reference
    std::vector<std::uint8_t> surface;      // SurfaceId per pixel
};

SyntheticDepth render_synthetic_depth(const CameraIntrinsics& K, const SyntheticSceneOptions& opt = {});

} // namespace barnes_hut::pc
