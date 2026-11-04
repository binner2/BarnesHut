#include "point_cloud.h"
#include "../stdinc.h"

#include <algorithm>
#include <bit>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <sstream>

#ifdef BH_HAVE_PNG
#include <png.h>
#endif

namespace barnes_hut::pc {

namespace {

std::string lower_extension(std::string_view filename) {
    const auto dot = filename.find_last_of('.');
    if (dot == std::string_view::npos) return {};
    std::string ext(filename.substr(dot + 1));
    for (auto& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext;
}

// Skip whitespace and '#' comments in a PNM header.
bool next_header_token(std::istream& is, std::string& token) {
    token.clear();
    int c;
    while ((c = is.peek()) != EOF) {
        if (c == '#') {
            std::string line;
            std::getline(is, line);
        } else if (std::isspace(c)) {
            is.get();
        } else {
            break;
        }
    }
    return static_cast<bool>(is >> token);
}

constexpr bool host_is_little_endian() { return std::endian::native == std::endian::little; }

template <class T>
T byteswap_value(T v) {
    std::array<std::uint8_t, sizeof(T)> bytes;
    std::memcpy(bytes.data(), &v, sizeof(T));
    std::reverse(bytes.begin(), bytes.end());
    std::memcpy(&v, bytes.data(), sizeof(T));
    return v;
}

} // namespace

// ---------------------------------------------------------------------------
// PointCloud
// ---------------------------------------------------------------------------

void PointCloud::reserve(Index n) {
    points.reserve(n);
}

void PointCloud::clear() {
    points.clear();
    normals.clear();
    colors.clear();
    scalars.clear();
}

PointCloud PointCloud::select(std::span<const Index> indices) const {
    PointCloud out;
    out.points.reserve(indices.size());
    const bool n = has_normals(), c = has_colors(), s = has_scalars();
    if (n) out.normals.reserve(indices.size());
    if (c) out.colors.reserve(indices.size());
    if (s) out.scalars.reserve(indices.size());
    for (const Index i : indices) {
        out.points.push_back(points[i]);
        if (n) out.normals.push_back(normals[i]);
        if (c) out.colors.push_back(colors[i]);
        if (s) out.scalars.push_back(scalars[i]);
    }
    return out;
}

bool PointCloud::bounding_box(Vector3D& lo, Vector3D& hi) const noexcept {
    if (points.empty()) return false;
    lo = hi = points.front();
    for (const auto& p : points)
        for (int k = 0; k < 3; ++k) {
            lo[k] = std::min(lo[k], p[k]);
            hi[k] = std::max(hi[k], p[k]);
        }
    return true;
}

Vector3D PointCloud::centroid() const noexcept {
    Vector3D acc{0.0};
    if (points.empty()) return acc;
    for (const auto& p : points) acc += p;
    return acc / static_cast<Real>(points.size());
}

// ---------------------------------------------------------------------------
// Depth -> cloud
// ---------------------------------------------------------------------------

PointCloud depth_to_point_cloud(const DepthImage& depth, const CameraIntrinsics& K,
                                const DepthToCloudOptions& opt,
                                std::vector<std::pair<int, int>>* pixel_of_point) {
    PointCloud cloud;
    const int stride = std::max(1, opt.stride);
    cloud.reserve(static_cast<Index>(depth.width) * depth.height / (stride * stride));
    if (pixel_of_point) pixel_of_point->clear();
    for (int v = 0; v < depth.height; v += stride) {
        for (int u = 0; u < depth.width; u += stride) {
            const std::uint16_t raw = depth.at(u, v);
            if (raw == 0) continue;
            const Real z = static_cast<Real>(raw) / K.depth_scale;
            if (z < opt.min_depth || z > opt.max_depth) continue;
            cloud.points.push_back(K.back_project(static_cast<Real>(u), static_cast<Real>(v), z));
            if (pixel_of_point) pixel_of_point->emplace_back(u, v);
        }
    }
    return cloud;
}

// ---------------------------------------------------------------------------
// Depth readers
// ---------------------------------------------------------------------------

std::optional<DepthImage> load_depth_pgm(std::string_view filename) {
    std::ifstream in(std::string(filename), std::ios::binary);
    if (!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return std::nullopt;
    }
    std::string magic, tok;
    if (!next_header_token(in, magic) || (magic != "P2" && magic != "P5")) {
        std::cerr << "Error: " << filename << " is not a PGM (P2/P5) file\n";
        return std::nullopt;
    }
    DepthImage img;
    int maxval = 0;
    if (!next_header_token(in, tok)) return std::nullopt;
    img.width = std::stoi(tok);
    if (!next_header_token(in, tok)) return std::nullopt;
    img.height = std::stoi(tok);
    if (!next_header_token(in, tok)) return std::nullopt;
    maxval = std::stoi(tok);
    if (img.width <= 0 || img.height <= 0 || maxval <= 0 || maxval > 65535) return std::nullopt;
    const std::size_t n = static_cast<std::size_t>(img.width) * img.height;
    img.data.resize(n);

    if (magic == "P2") {
        for (std::size_t i = 0; i < n; ++i) {
            long v;
            if (!(in >> v)) return std::nullopt;
            img.data[i] = static_cast<std::uint16_t>(std::clamp<long>(v, 0, 65535));
        }
        return img;
    }
    in.get();   // single whitespace after maxval
    if (maxval < 256) {
        std::vector<std::uint8_t> buf(n);
        in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(n));
        if (!in) return std::nullopt;
        for (std::size_t i = 0; i < n; ++i) img.data[i] = buf[i];
    } else {
        std::vector<std::uint8_t> buf(2 * n);
        in.read(reinterpret_cast<char*>(buf.data()), static_cast<std::streamsize>(2 * n));
        if (!in) return std::nullopt;
        for (std::size_t i = 0; i < n; ++i)   // PGM 16-bit samples are big-endian
            img.data[i] = static_cast<std::uint16_t>((buf[2 * i] << 8) | buf[2 * i + 1]);
    }
    return img;
}

std::optional<DepthImage> load_depth_raw(std::string_view filename, int width, int height) {
    std::ifstream in(std::string(filename), std::ios::binary | std::ios::ate);
    if (!in || width <= 0 || height <= 0) return std::nullopt;
    const std::size_t n = static_cast<std::size_t>(width) * height;
    if (static_cast<std::size_t>(in.tellg()) != 2 * n) {
        std::cerr << "Error: raw depth size mismatch for " << filename << "\n";
        return std::nullopt;
    }
    in.seekg(0);
    DepthImage img;
    img.width = width;
    img.height = height;
    img.data.resize(n);
    in.read(reinterpret_cast<char*>(img.data.data()), static_cast<std::streamsize>(2 * n));
    if (!host_is_little_endian())
        for (auto& v : img.data) v = byteswap_value(v);
    return img;
}

bool png_support_compiled() noexcept {
#ifdef BH_HAVE_PNG
    return true;
#else
    return false;
#endif
}

std::optional<DepthImage> load_depth_png(std::string_view filename) {
#ifndef BH_HAVE_PNG
    std::cerr << "Error: PNG support not compiled in (build with libpng, -DBH_HAVE_PNG). "
                 "Convert with tools/depth_png_to_pgm.py or use PGM/raw input: " << filename << "\n";
    return std::nullopt;
#else
    FILE* fp = std::fopen(std::string(filename).c_str(), "rb");
    if (!fp) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return std::nullopt;
    }
    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop info = png ? png_create_info_struct(png) : nullptr;
    if (!png || !info) {
        if (png) png_destroy_read_struct(&png, info ? &info : nullptr, nullptr);
        std::fclose(fp);
        return std::nullopt;
    }
    if (setjmp(png_jmpbuf(png))) {
        png_destroy_read_struct(&png, &info, nullptr);
        std::fclose(fp);
        std::cerr << "Error: libpng failed on " << filename << "\n";
        return std::nullopt;
    }
    png_init_io(png, fp);
    png_read_info(png, info);
    const png_uint_32 width = png_get_image_width(png, info);
    const png_uint_32 height = png_get_image_height(png, info);
    const png_byte color_type = png_get_color_type(png, info);
    const png_byte bit_depth = png_get_bit_depth(png, info);

    if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
    if (color_type & PNG_COLOR_MASK_ALPHA) png_set_strip_alpha(png);
    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_RGB_ALPHA ||
        color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_rgb_to_gray_fixed(png, 1, -1, -1);
    if (bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
    if (bit_depth == 16 && host_is_little_endian()) png_set_swap(png);
    png_read_update_info(png, info);

    const png_size_t rowbytes = png_get_rowbytes(png, info);
    std::vector<png_byte> buffer(rowbytes * height);
    std::vector<png_bytep> rows(height);
    for (png_uint_32 y = 0; y < height; ++y) rows[y] = buffer.data() + y * rowbytes;
    png_read_image(png, rows.data());
    png_destroy_read_struct(&png, &info, nullptr);
    std::fclose(fp);

    DepthImage img;
    img.width = static_cast<int>(width);
    img.height = static_cast<int>(height);
    img.data.resize(static_cast<std::size_t>(width) * height);
    if (rowbytes == 2 * width) {
        std::memcpy(img.data.data(), buffer.data(), buffer.size());
    } else if (rowbytes == width) {
        for (std::size_t i = 0; i < img.data.size(); ++i) img.data[i] = buffer[i];
    } else {
        std::cerr << "Error: unsupported PNG layout in " << filename << " (rowbytes=" << rowbytes << ")\n";
        return std::nullopt;
    }
    return img;
#endif
}

std::optional<DepthImage> load_depth_image(std::string_view filename) {
    const std::string ext = lower_extension(filename);
    if (ext == "pgm") return load_depth_pgm(filename);
    if (ext == "png") return load_depth_png(filename);
    std::cerr << "Error: unknown depth image extension '" << ext << "' (use .png or .pgm)\n";
    return std::nullopt;
}

bool save_depth_pgm(const DepthImage& depth, std::string_view filename) {
    std::ofstream out(std::string(filename), std::ios::binary);
    if (!out) return false;
    out << "P5\n" << depth.width << " " << depth.height << "\n65535\n";
    std::vector<std::uint8_t> buf(depth.data.size() * 2);
    for (std::size_t i = 0; i < depth.data.size(); ++i) {
        buf[2 * i] = static_cast<std::uint8_t>(depth.data[i] >> 8);
        buf[2 * i + 1] = static_cast<std::uint8_t>(depth.data[i] & 0xFF);
    }
    out.write(reinterpret_cast<const char*>(buf.data()), static_cast<std::streamsize>(buf.size()));
    return static_cast<bool>(out);
}

std::optional<CameraIntrinsics> load_intrinsics(std::string_view filename) {
    std::ifstream in{std::string(filename)};
    if (!in) {
        std::cerr << "Error: cannot open intrinsics " << filename << "\n";
        return std::nullopt;
    }
    std::vector<Real> vals;
    std::string line;
    while (std::getline(in, line)) {
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.erase(hash);
        std::istringstream ls(line);
        Real v;
        while (ls >> v) vals.push_back(v);
    }
    if (vals.size() < 4) {
        std::cerr << "Error: intrinsics need at least fx fy cx cy\n";
        return std::nullopt;
    }
    CameraIntrinsics K;
    K.fx = vals[0];
    K.fy = vals[1];
    K.cx = vals[2];
    K.cy = vals[3];
    if (vals.size() >= 6) {
        K.width = static_cast<int>(vals[4]);
        K.height = static_cast<int>(vals[5]);
    }
    if (vals.size() >= 7) K.depth_scale = vals[6];
    return K;
}

// ---------------------------------------------------------------------------
// XYZ
// ---------------------------------------------------------------------------

std::optional<PointCloud> load_xyz(std::string_view filename) {
    std::ifstream in{std::string(filename)};
    if (!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return std::nullopt;
    }
    PointCloud cloud;
    std::string line;
    std::vector<Real> vals;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '/') continue;
        for (auto& c : line) if (c == ',' || c == ';') c = ' ';
        std::istringstream ls(line);
        vals.clear();
        Real v;
        while (ls >> v) vals.push_back(v);
        if (vals.size() < 3) continue;
        cloud.points.emplace_back(vals[0], vals[1], vals[2]);
        if (vals.size() >= 6) cloud.normals.emplace_back(vals[3], vals[4], vals[5]);
        if (vals.size() >= 9)
            cloud.colors.push_back(Color{static_cast<std::uint8_t>(std::clamp(vals[6], 0.0, 255.0)),
                                         static_cast<std::uint8_t>(std::clamp(vals[7], 0.0, 255.0)),
                                         static_cast<std::uint8_t>(std::clamp(vals[8], 0.0, 255.0))});
    }
    if (cloud.normals.size() != cloud.points.size()) cloud.normals.clear();
    if (cloud.colors.size() != cloud.points.size()) cloud.colors.clear();
    return cloud;
}

bool save_xyz(const PointCloud& cloud, std::string_view filename) {
    std::ofstream out{std::string(filename)};
    if (!out) return false;
    out << std::setprecision(9);
    for (Index i = 0; i < cloud.size(); ++i) {
        const auto& p = cloud.points[i];
        out << p.x() << ' ' << p.y() << ' ' << p.z();
        if (cloud.has_normals()) {
            const auto& n = cloud.normals[i];
            out << ' ' << n.x() << ' ' << n.y() << ' ' << n.z();
        }
        if (cloud.has_colors()) {
            const auto& c = cloud.colors[i];
            out << ' ' << int(c.r) << ' ' << int(c.g) << ' ' << int(c.b);
        }
        if (cloud.has_scalars()) out << ' ' << cloud.scalars[i];
        out << '\n';
    }
    return static_cast<bool>(out);
}

// ---------------------------------------------------------------------------
// PLY
// ---------------------------------------------------------------------------

namespace {

enum class PlyType { Int8, UInt8, Int16, UInt16, Int32, UInt32, Float32, Float64, Unknown };

PlyType ply_type_from(const std::string& s) {
    if (s == "char" || s == "int8") return PlyType::Int8;
    if (s == "uchar" || s == "uint8") return PlyType::UInt8;
    if (s == "short" || s == "int16") return PlyType::Int16;
    if (s == "ushort" || s == "uint16") return PlyType::UInt16;
    if (s == "int" || s == "int32") return PlyType::Int32;
    if (s == "uint" || s == "uint32") return PlyType::UInt32;
    if (s == "float" || s == "float32") return PlyType::Float32;
    if (s == "double" || s == "float64") return PlyType::Float64;
    return PlyType::Unknown;
}

struct PlyProperty {
    std::string name;
    PlyType type = PlyType::Unknown;
    bool is_list = false;
    PlyType count_type = PlyType::Unknown;
};

struct PlyElement {
    std::string name;
    std::size_t count = 0;
    std::vector<PlyProperty> properties;
};

template <class T>
double read_binary_as_double(std::istream& in, bool swap) {
    T v{};
    in.read(reinterpret_cast<char*>(&v), sizeof(T));
    if (swap) v = byteswap_value(v);
    return static_cast<double>(v);
}

double read_ply_scalar(std::istream& in, PlyType t, bool binary, bool swap) {
    if (!binary) {
        double v = 0.0;
        in >> v;
        return v;
    }
    switch (t) {
        case PlyType::Int8: return read_binary_as_double<std::int8_t>(in, swap);
        case PlyType::UInt8: return read_binary_as_double<std::uint8_t>(in, swap);
        case PlyType::Int16: return read_binary_as_double<std::int16_t>(in, swap);
        case PlyType::UInt16: return read_binary_as_double<std::uint16_t>(in, swap);
        case PlyType::Int32: return read_binary_as_double<std::int32_t>(in, swap);
        case PlyType::UInt32: return read_binary_as_double<std::uint32_t>(in, swap);
        case PlyType::Float32: return read_binary_as_double<float>(in, swap);
        case PlyType::Float64: return read_binary_as_double<double>(in, swap);
        default: return 0.0;
    }
}

} // namespace

std::optional<PointCloud> load_ply(std::string_view filename) {
    std::ifstream in(std::string(filename), std::ios::binary);
    if (!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return std::nullopt;
    }
    std::string line;
    std::getline(in, line);
    if (line.rfind("ply", 0) != 0) {
        std::cerr << "Error: " << filename << " is not a PLY file\n";
        return std::nullopt;
    }
    bool binary = false, swap = false;
    std::vector<PlyElement> elements;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        std::istringstream ls(line);
        std::string key;
        ls >> key;
        if (key == "format") {
            std::string fmt;
            ls >> fmt;
            if (fmt == "ascii") binary = false;
            else if (fmt == "binary_little_endian") { binary = true; swap = !host_is_little_endian(); }
            else if (fmt == "binary_big_endian") { binary = true; swap = host_is_little_endian(); }
            else return std::nullopt;
        } else if (key == "element") {
            PlyElement e;
            ls >> e.name >> e.count;
            elements.push_back(e);
        } else if (key == "property") {
            if (elements.empty()) return std::nullopt;
            PlyProperty p;
            std::string t;
            ls >> t;
            if (t == "list") {
                std::string ct, vt;
                ls >> ct >> vt >> p.name;
                p.is_list = true;
                p.count_type = ply_type_from(ct);
                p.type = ply_type_from(vt);
            } else {
                p.type = ply_type_from(t);
                ls >> p.name;
            }
            elements.back().properties.push_back(p);
        } else if (key == "end_header") {
            break;
        }
    }

    PointCloud cloud;
    for (const PlyElement& e : elements) {
        if (e.name != "vertex") {
            // Skip foreign elements (needed only when they precede 'vertex').
            for (std::size_t i = 0; i < e.count; ++i)
                for (const PlyProperty& p : e.properties) {
                    if (p.is_list) {
                        const auto n = static_cast<std::size_t>(read_ply_scalar(in, p.count_type, binary, swap));
                        for (std::size_t j = 0; j < n; ++j) read_ply_scalar(in, p.type, binary, swap);
                    } else {
                        read_ply_scalar(in, p.type, binary, swap);
                    }
                }
            continue;
        }
        int ix = -1, iy = -1, iz = -1, inx = -1, iny = -1, inz = -1, ir = -1, ig = -1, ib = -1, is = -1;
        for (int k = 0; k < static_cast<int>(e.properties.size()); ++k) {
            const std::string& n = e.properties[k].name;
            if (n == "x") ix = k; else if (n == "y") iy = k; else if (n == "z") iz = k;
            else if (n == "nx") inx = k; else if (n == "ny") iny = k; else if (n == "nz") inz = k;
            else if (n == "red" || n == "r" || n == "diffuse_red") ir = k;
            else if (n == "green" || n == "g" || n == "diffuse_green") ig = k;
            else if (n == "blue" || n == "b" || n == "diffuse_blue") ib = k;
            else if (is < 0 && (n == "scalar" || n == "intensity" || n == "density" || n == "label" || n == "value" ||
                                n.rfind("scalar_", 0) == 0)) is = k;
        }
        if (ix < 0 || iy < 0 || iz < 0) {
            std::cerr << "Error: PLY vertex element lacks x/y/z\n";
            return std::nullopt;
        }
        const bool has_n = inx >= 0 && iny >= 0 && inz >= 0;
        const bool has_c = ir >= 0 && ig >= 0 && ib >= 0;
        cloud.points.reserve(e.count);
        std::vector<double> vals(e.properties.size());
        for (std::size_t i = 0; i < e.count; ++i) {
            for (std::size_t k = 0; k < e.properties.size(); ++k) {
                const PlyProperty& p = e.properties[k];
                if (p.is_list) {
                    const auto n = static_cast<std::size_t>(read_ply_scalar(in, p.count_type, binary, swap));
                    for (std::size_t j = 0; j < n; ++j) read_ply_scalar(in, p.type, binary, swap);
                    vals[k] = 0.0;
                } else {
                    vals[k] = read_ply_scalar(in, p.type, binary, swap);
                }
            }
            if (!in) {
                std::cerr << "Error: truncated PLY after " << i << " vertices\n";
                return std::nullopt;
            }
            cloud.points.emplace_back(vals[ix], vals[iy], vals[iz]);
            if (has_n) cloud.normals.emplace_back(vals[inx], vals[iny], vals[inz]);
            if (has_c) {
                auto to8 = [&](int k) {
                    const PlyType t = e.properties[k].type;
                    double v = vals[k];
                    if (t == PlyType::Float32 || t == PlyType::Float64) v *= 255.0;
                    return static_cast<std::uint8_t>(std::clamp(v, 0.0, 255.0));
                };
                cloud.colors.push_back(Color{to8(ir), to8(ig), to8(ib)});
            }
            if (is >= 0) cloud.scalars.push_back(vals[is]);
        }
        break;  // elements after vertex are not needed
    }
    if (cloud.points.empty()) {
        std::cerr << "Error: PLY has no vertex element\n";
        return std::nullopt;
    }
    return cloud;
}

bool save_ply(const PointCloud& cloud, std::string_view filename, bool binary) {
    std::ofstream out(std::string(filename), std::ios::binary);
    if (!out) return false;
    const bool n = cloud.has_normals(), c = cloud.has_colors(), s = cloud.has_scalars();
    out << "ply\nformat " << (binary ? "binary_little_endian" : "ascii") << " 1.0\n"
        << "comment generated by BarnesHut/pointcloud\n"
        << "element vertex " << cloud.size() << "\n"
        << "property float x\nproperty float y\nproperty float z\n";
    if (n) out << "property float nx\nproperty float ny\nproperty float nz\n";
    if (c) out << "property uchar red\nproperty uchar green\nproperty uchar blue\n";
    if (s) out << "property float scalar\n";
    out << "end_header\n";

    if (binary) {
        std::vector<char> buf;
        buf.reserve(cloud.size() * 32);
        auto put_f = [&](Real v) {
            float f = static_cast<float>(v);
            if (!host_is_little_endian()) f = byteswap_value(f);
            const char* b = reinterpret_cast<const char*>(&f);
            buf.insert(buf.end(), b, b + 4);
        };
        for (Index i = 0; i < cloud.size(); ++i) {
            put_f(cloud.points[i].x()); put_f(cloud.points[i].y()); put_f(cloud.points[i].z());
            if (n) { put_f(cloud.normals[i].x()); put_f(cloud.normals[i].y()); put_f(cloud.normals[i].z()); }
            if (c) { buf.push_back(static_cast<char>(cloud.colors[i].r)); buf.push_back(static_cast<char>(cloud.colors[i].g)); buf.push_back(static_cast<char>(cloud.colors[i].b)); }
            if (s) put_f(cloud.scalars[i]);
        }
        out.write(buf.data(), static_cast<std::streamsize>(buf.size()));
    } else {
        out << std::setprecision(8);
        for (Index i = 0; i < cloud.size(); ++i) {
            out << cloud.points[i].x() << ' ' << cloud.points[i].y() << ' ' << cloud.points[i].z();
            if (n) out << ' ' << cloud.normals[i].x() << ' ' << cloud.normals[i].y() << ' ' << cloud.normals[i].z();
            if (c) out << ' ' << int(cloud.colors[i].r) << ' ' << int(cloud.colors[i].g) << ' ' << int(cloud.colors[i].b);
            if (s) out << ' ' << cloud.scalars[i];
            out << '\n';
        }
    }
    return static_cast<bool>(out);
}

// ---------------------------------------------------------------------------
// PCD
// ---------------------------------------------------------------------------

std::optional<PointCloud> load_pcd(std::string_view filename) {
    std::ifstream in(std::string(filename), std::ios::binary);
    if (!in) {
        std::cerr << "Error: cannot open " << filename << "\n";
        return std::nullopt;
    }
    std::vector<std::string> fields, types;
    std::vector<int> sizes, counts;
    std::size_t points = 0;
    std::string data_mode, line;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ls(line);
        std::string key;
        ls >> key;
        if (key == "FIELDS") { std::string f; while (ls >> f) fields.push_back(f); }
        else if (key == "SIZE") { int v; while (ls >> v) sizes.push_back(v); }
        else if (key == "TYPE") { std::string t; while (ls >> t) types.push_back(t); }
        else if (key == "COUNT") { int v; while (ls >> v) counts.push_back(v); }
        else if (key == "POINTS") { ls >> points; }
        else if (key == "WIDTH" && points == 0) { ls >> points; }
        else if (key == "HEIGHT" && points > 0) { std::size_t h; if (ls >> h) { /* WIDTH already read */ } }
        else if (key == "DATA") { ls >> data_mode; break; }
    }
    if (fields.empty() || sizes.size() != fields.size() || types.size() != fields.size()) {
        std::cerr << "Error: malformed PCD header in " << filename << "\n";
        return std::nullopt;
    }
    if (counts.size() != fields.size()) counts.assign(fields.size(), 1);
    if (data_mode == "binary_compressed") {
        std::cerr << "Error: binary_compressed PCD is not supported; re-save as binary or ascii\n";
        return std::nullopt;
    }
    const bool binary = data_mode == "binary";

    // Field offsets (binary) / token indices (ascii).
    std::vector<std::size_t> offsets(fields.size());
    std::size_t stride = 0;
    std::vector<std::size_t> token_index(fields.size());
    std::size_t tokens = 0;
    for (std::size_t k = 0; k < fields.size(); ++k) {
        offsets[k] = stride;
        stride += static_cast<std::size_t>(sizes[k]) * counts[k];
        token_index[k] = tokens;
        tokens += counts[k];
    }
    auto find = [&](std::string_view name) -> int {
        for (std::size_t k = 0; k < fields.size(); ++k)
            if (fields[k] == name) return static_cast<int>(k);
        return -1;
    };
    const int ix = find("x"), iy = find("y"), iz = find("z");
    const int inx = find("normal_x"), iny = find("normal_y"), inz = find("normal_z");
    const int irgb = find("rgb") >= 0 ? find("rgb") : find("rgba");
    int is = find("intensity");
    if (is < 0) is = find("label");
    if (ix < 0 || iy < 0 || iz < 0) {
        std::cerr << "Error: PCD lacks x/y/z fields\n";
        return std::nullopt;
    }
    const bool has_n = inx >= 0 && iny >= 0 && inz >= 0;

    PointCloud cloud;
    cloud.points.reserve(points);

    auto decode_rgb = [](float f) {
        std::uint32_t u;
        std::memcpy(&u, &f, 4);
        return Color{static_cast<std::uint8_t>((u >> 16) & 0xFF), static_cast<std::uint8_t>((u >> 8) & 0xFF),
                     static_cast<std::uint8_t>(u & 0xFF)};
    };

    if (binary) {
        std::vector<char> rec(stride);
        auto get = [&](int k) -> double {
            const char* p = rec.data() + offsets[k];
            const int sz = sizes[k];
            const char t = types[k][0];
            if (t == 'F' && sz == 4) { float v; std::memcpy(&v, p, 4); return v; }
            if (t == 'F' && sz == 8) { double v; std::memcpy(&v, p, 8); return v; }
            if (t == 'U' && sz == 1) { std::uint8_t v; std::memcpy(&v, p, 1); return v; }
            if (t == 'U' && sz == 2) { std::uint16_t v; std::memcpy(&v, p, 2); return v; }
            if (t == 'U' && sz == 4) { std::uint32_t v; std::memcpy(&v, p, 4); return v; }
            if (t == 'I' && sz == 1) { std::int8_t v; std::memcpy(&v, p, 1); return v; }
            if (t == 'I' && sz == 2) { std::int16_t v; std::memcpy(&v, p, 2); return v; }
            if (t == 'I' && sz == 4) { std::int32_t v; std::memcpy(&v, p, 4); return v; }
            return 0.0;
        };
        for (std::size_t i = 0; i < points; ++i) {
            in.read(rec.data(), static_cast<std::streamsize>(stride));
            if (!in) break;
            const Vector3D p{get(ix), get(iy), get(iz)};
            if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z())) continue;
            cloud.points.push_back(p);
            if (has_n) cloud.normals.emplace_back(get(inx), get(iny), get(inz));
            if (irgb >= 0) {
                if (types[irgb][0] == 'F') { float f; std::memcpy(&f, rec.data() + offsets[irgb], 4); cloud.colors.push_back(decode_rgb(f)); }
                else { const auto u = static_cast<std::uint32_t>(get(irgb)); cloud.colors.push_back(Color{static_cast<std::uint8_t>((u >> 16) & 0xFF), static_cast<std::uint8_t>((u >> 8) & 0xFF), static_cast<std::uint8_t>(u & 0xFF)}); }
            }
            if (is >= 0) cloud.scalars.push_back(get(is));
        }
    } else {
        std::vector<std::string> tok;
        while (std::getline(in, line) && cloud.points.size() < points) {
            std::istringstream ls(line);
            tok.assign(std::istream_iterator<std::string>(ls), std::istream_iterator<std::string>());
            if (tok.size() < tokens) continue;
            auto get = [&](int k) { return std::stod(tok[token_index[k]]); };
            const Vector3D p{get(ix), get(iy), get(iz)};
            if (!std::isfinite(p.x()) || !std::isfinite(p.y()) || !std::isfinite(p.z())) continue;
            cloud.points.push_back(p);
            if (has_n) cloud.normals.emplace_back(get(inx), get(iny), get(inz));
            if (irgb >= 0) {
                if (types[irgb][0] == 'F') cloud.colors.push_back(decode_rgb(static_cast<float>(std::stod(tok[token_index[irgb]]))));
                else { const auto u = static_cast<std::uint32_t>(std::stoul(tok[token_index[irgb]])); cloud.colors.push_back(Color{static_cast<std::uint8_t>((u >> 16) & 0xFF), static_cast<std::uint8_t>((u >> 8) & 0xFF), static_cast<std::uint8_t>(u & 0xFF)}); }
            }
            if (is >= 0) cloud.scalars.push_back(get(is));
        }
    }
    if (cloud.normals.size() != cloud.points.size()) cloud.normals.clear();
    if (cloud.colors.size() != cloud.points.size()) cloud.colors.clear();
    if (cloud.scalars.size() != cloud.points.size()) cloud.scalars.clear();
    return cloud;
}

// ---------------------------------------------------------------------------
// Dispatch
// ---------------------------------------------------------------------------

std::optional<PointCloud> load_point_cloud(std::string_view filename) {
    const std::string ext = lower_extension(filename);
    if (ext == "ply") return load_ply(filename);
    if (ext == "pcd") return load_pcd(filename);
    if (ext == "xyz" || ext == "txt" || ext == "csv" || ext == "asc" || ext == "pts") return load_xyz(filename);
    std::cerr << "Error: unknown point cloud extension '" << ext << "' (ply, pcd, xyz, txt)\n";
    return std::nullopt;
}

bool save_point_cloud(const PointCloud& cloud, std::string_view filename) {
    const std::string ext = lower_extension(filename);
    if (ext == "ply") return save_ply(cloud, filename, true);
    if (ext == "xyz" || ext == "txt" || ext == "csv" || ext == "asc") return save_xyz(cloud, filename);
    std::cerr << "Error: unsupported output extension '" << ext << "' (ply, xyz, txt)\n";
    return false;
}

// ---------------------------------------------------------------------------
// Colour helpers
// ---------------------------------------------------------------------------

void colorize_by_scalar(PointCloud& cloud) {
    if (!cloud.has_scalars()) return;
    // Sampled viridis anchors (Matplotlib), interpolated linearly.
    static constexpr std::array<std::array<double, 3>, 9> anchors{{
        {68, 1, 84}, {72, 36, 117}, {65, 68, 135}, {53, 95, 141}, {42, 120, 142},
        {33, 145, 140}, {53, 183, 121}, {144, 215, 67}, {253, 231, 37}}};
    std::vector<Real> sorted = cloud.scalars;
    std::sort(sorted.begin(), sorted.end());
    const Real lo = sorted[static_cast<std::size_t>(0.02 * (sorted.size() - 1))];
    const Real hi = sorted[static_cast<std::size_t>(0.98 * (sorted.size() - 1))];
    const Real span = hi > lo ? hi - lo : 1.0;
    cloud.colors.resize(cloud.size());
    for (Index i = 0; i < cloud.size(); ++i) {
        const double t = std::clamp((cloud.scalars[i] - lo) / span, 0.0, 1.0) * (anchors.size() - 1);
        const std::size_t a = std::min<std::size_t>(static_cast<std::size_t>(t), anchors.size() - 2);
        const double f = t - static_cast<double>(a);
        Color c;
        c.r = static_cast<std::uint8_t>(anchors[a][0] + f * (anchors[a + 1][0] - anchors[a][0]));
        c.g = static_cast<std::uint8_t>(anchors[a][1] + f * (anchors[a + 1][1] - anchors[a][1]));
        c.b = static_cast<std::uint8_t>(anchors[a][2] + f * (anchors[a + 1][2] - anchors[a][2]));
        cloud.colors[i] = c;
    }
}

void colorize_by_label(PointCloud& cloud) {
    if (!cloud.has_scalars()) return;
    cloud.colors.resize(cloud.size());
    for (Index i = 0; i < cloud.size(); ++i) {
        const long label = static_cast<long>(std::llround(cloud.scalars[i]));
        if (label < 0) {
            cloud.colors[i] = Color{128, 128, 128};
            continue;
        }
        const double hue = std::fmod(static_cast<double>(label) * 0.6180339887498949, 1.0) * 6.0;
        const double s = 0.65, v = 0.95;
        const int sector = static_cast<int>(hue);
        const double f = hue - sector;
        const double p = v * (1 - s), q = v * (1 - s * f), t = v * (1 - s * (1 - f));
        double r, g, b;
        switch (sector % 6) {
            case 0: r = v; g = t; b = p; break;
            case 1: r = q; g = v; b = p; break;
            case 2: r = p; g = v; b = t; break;
            case 3: r = p; g = q; b = v; break;
            case 4: r = t; g = p; b = v; break;
            default: r = v; g = p; b = q; break;
        }
        cloud.colors[i] = Color{static_cast<std::uint8_t>(r * 255), static_cast<std::uint8_t>(g * 255),
                                static_cast<std::uint8_t>(b * 255)};
    }
}

} // namespace barnes_hut::pc
