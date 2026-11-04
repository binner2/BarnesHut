#include "synthetic.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace barnes_hut::pc {

namespace {

struct Hit {
    Real t = std::numeric_limits<Real>::infinity();
    SurfaceId id = SurfaceId::None;
};

void consider(Hit& best, Real t, SurfaceId id) {
    if (t > 0.0 && t < best.t) {
        best.t = t;
        best.id = id;
    }
}

} // namespace

SyntheticDepth render_synthetic_depth(const CameraIntrinsics& K, const SyntheticSceneOptions& opt) {
    SyntheticDepth out;
    out.depth.width = out.depth_clean.width = K.width;
    out.depth.height = out.depth_clean.height = K.height;
    const std::size_t n = static_cast<std::size_t>(K.width) * K.height;
    out.depth.data.assign(n, 0);
    out.depth_clean.data.assign(n, 0);
    out.surface.assign(n, 0);

    std::mt19937_64 rng(opt.seed);
    std::normal_distribution<Real> gauss(0.0, 1.0);
    std::uniform_real_distribution<Real> uni(0.0, 1.0);

    for (int v = 0; v < K.height; ++v) {
        for (int u = 0; u < K.width; ++u) {
            // Ray through pixel centre, parametrised by z (dir has z = 1).
            const Vector3D dir{(u - K.cx) / K.fx, (v - K.cy) / K.fy, 1.0};
            Hit hit;

            // Back wall z = wall_distance.
            consider(hit, opt.wall_distance, SurfaceId::BackWall);
            // Floor y = +floor_height (image y axis points down).
            if (dir.y() > 1e-9) consider(hit, opt.floor_height / dir.y(), SurfaceId::Floor);
            // Side wall x = +side_wall_x.
            if (dir.x() > 1e-9) consider(hit, opt.side_wall_x / dir.x(), SurfaceId::SideWall);
            // Sphere.
            {
                const Vector3D oc = -opt.sphere_center;          // origin - centre
                const Real a = dir.dot(dir), b = 2.0 * oc.dot(dir);
                const Real c = oc.dot(oc) - opt.sphere_radius * opt.sphere_radius;
                const Real disc = b * b - 4 * a * c;
                if (disc >= 0.0) consider(hit, (-b - std::sqrt(disc)) / (2 * a), SurfaceId::Sphere);
            }
            // Axis-aligned box (slab test).
            {
                Real t_near = 0.0, t_far = std::numeric_limits<Real>::infinity();
                bool ok = true;
                for (int k = 0; k < 3 && ok; ++k) {
                    if (std::abs(dir[k]) < 1e-12) {
                        if (0.0 < opt.box_min[k] || 0.0 > opt.box_max[k]) ok = false;
                        continue;
                    }
                    Real t0 = opt.box_min[k] / dir[k], t1 = opt.box_max[k] / dir[k];
                    if (t0 > t1) std::swap(t0, t1);
                    t_near = std::max(t_near, t0);
                    t_far = std::min(t_far, t1);
                    if (t_near > t_far) ok = false;
                }
                if (ok) consider(hit, t_near, SurfaceId::Box);
            }

            if (hit.id == SurfaceId::None) continue;
            const Real z = hit.t;   // dir.z == 1, so t is the depth
            const std::size_t idx = static_cast<std::size_t>(v) * K.width + u;
            out.surface[idx] = static_cast<std::uint8_t>(hit.id);
            const Real clean_raw = std::round(z * K.depth_scale);
            out.depth_clean.data[idx] = static_cast<std::uint16_t>(std::clamp<Real>(clean_raw, 0.0, 65535.0));
            if (uni(rng) < opt.dropout_fraction) continue;
            const Real z_noisy = z + opt.noise_coefficient * z * z * gauss(rng);
            const Real raw = std::round(z_noisy * K.depth_scale);
            out.depth.data[idx] = static_cast<std::uint16_t>(std::clamp<Real>(raw, 0.0, 65535.0));
        }
    }
    return out;
}

} // namespace barnes_hut::pc
