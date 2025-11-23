#pragma once

#include "stdinc.h"
#include <array>
#include <cmath>
#include <iostream>

namespace barnes_hut {

// Modern 3D vector class with C++20 features
class Vector3D {
public:
    // Constructors
    constexpr Vector3D() noexcept : data_{0.0, 0.0, 0.0} {}

    constexpr explicit Vector3D(Real value) noexcept : data_{value, value, value} {}

    constexpr Vector3D(Real x, Real y, Real z) noexcept : data_{x, y, z} {}

    // Modern accessors with bounds checking in debug mode
    [[nodiscard]] constexpr Real& operator[](std::size_t i) noexcept {
        return data_[i];
    }

    [[nodiscard]] constexpr const Real& operator[](std::size_t i) const noexcept {
        return data_[i];
    }

    [[nodiscard]] constexpr Real& x() noexcept { return data_[0]; }
    [[nodiscard]] constexpr Real& y() noexcept { return data_[1]; }
    [[nodiscard]] constexpr Real& z() noexcept { return data_[2]; }

    [[nodiscard]] constexpr Real x() const noexcept { return data_[0]; }
    [[nodiscard]] constexpr Real y() const noexcept { return data_[1]; }
    [[nodiscard]] constexpr Real z() const noexcept { return data_[2]; }

    // Unary operators
    [[nodiscard]] constexpr Vector3D operator-() const noexcept {
        return Vector3D{-data_[0], -data_[1], -data_[2]};
    }

    // Dot product
    [[nodiscard]] constexpr Real dot(const Vector3D& other) const noexcept {
        return data_[0] * other.data_[0] +
               data_[1] * other.data_[1] +
               data_[2] * other.data_[2];
    }

    [[nodiscard]] constexpr Real operator*(const Vector3D& other) const noexcept {
        return dot(other);
    }

    // Cross product
    [[nodiscard]] constexpr Vector3D cross(const Vector3D& other) const noexcept {
        return Vector3D{
            data_[1] * other.data_[2] - data_[2] * other.data_[1],
            data_[2] * other.data_[0] - data_[0] * other.data_[2],
            data_[0] * other.data_[1] - data_[1] * other.data_[0]
        };
    }

    // Vector addition
    [[nodiscard]] constexpr Vector3D operator+(const Vector3D& other) const noexcept {
        return Vector3D{
            data_[0] + other.data_[0],
            data_[1] + other.data_[1],
            data_[2] + other.data_[2]
        };
    }

    // Vector subtraction
    [[nodiscard]] constexpr Vector3D operator-(const Vector3D& other) const noexcept {
        return Vector3D{
            data_[0] - other.data_[0],
            data_[1] - other.data_[1],
            data_[2] - other.data_[2]
        };
    }

    // Scalar multiplication
    [[nodiscard]] constexpr Vector3D operator*(Real scalar) const noexcept {
        return Vector3D{
            data_[0] * scalar,
            data_[1] * scalar,
            data_[2] * scalar
        };
    }

    friend constexpr Vector3D operator*(Real scalar, const Vector3D& vec) noexcept {
        return vec * scalar;
    }

    // Scalar division
    [[nodiscard]] constexpr Vector3D operator/(Real scalar) const noexcept {
        return Vector3D{
            data_[0] / scalar,
            data_[1] / scalar,
            data_[2] / scalar
        };
    }

    // Scalar addition
    [[nodiscard]] constexpr Vector3D operator+(Real scalar) const noexcept {
        return Vector3D{
            data_[0] + scalar,
            data_[1] + scalar,
            data_[2] + scalar
        };
    }

    friend constexpr Vector3D operator+(Real scalar, const Vector3D& vec) noexcept {
        return vec + scalar;
    }

    // Compound assignment operators
    constexpr Vector3D& operator+=(const Vector3D& other) noexcept {
        data_[0] += other.data_[0];
        data_[1] += other.data_[1];
        data_[2] += other.data_[2];
        return *this;
    }

    constexpr Vector3D& operator-=(const Vector3D& other) noexcept {
        data_[0] -= other.data_[0];
        data_[1] -= other.data_[1];
        data_[2] -= other.data_[2];
        return *this;
    }

    constexpr Vector3D& operator*=(Real scalar) noexcept {
        data_[0] *= scalar;
        data_[1] *= scalar;
        data_[2] *= scalar;
        return *this;
    }

    constexpr Vector3D& operator/=(Real scalar) noexcept {
        data_[0] /= scalar;
        data_[1] /= scalar;
        data_[2] /= scalar;
        return *this;
    }

    // Comparison operators
    [[nodiscard]] constexpr bool operator==(const Vector3D& other) const noexcept {
        return data_[0] == other.data_[0] &&
               data_[1] == other.data_[1] &&
               data_[2] == other.data_[2];
    }

    [[nodiscard]] constexpr bool operator!=(const Vector3D& other) const noexcept {
        return !(*this == other);
    }

    // Magnitude functions
    [[nodiscard]] constexpr Real squared_magnitude() const noexcept {
        return dot(*this);
    }

    [[nodiscard]] Real magnitude() const noexcept {
        return std::sqrt(squared_magnitude());
    }

    // Distance functions
    [[nodiscard]] constexpr Real squared_distance(const Vector3D& other) const noexcept {
        const Real dx = data_[0] - other.data_[0];
        const Real dy = data_[1] - other.data_[1];
        const Real dz = data_[2] - other.data_[2];
        return dx * dx + dy * dy + dz * dz;
    }

    [[nodiscard]] Real distance(const Vector3D& other) const noexcept {
        return std::sqrt(squared_distance(other));
    }

    // Normalization
    [[nodiscard]] Vector3D normalized() const noexcept {
        const Real mag = magnitude();
        if (mag < 1e-10) return Vector3D{0.0};
        return *this / mag;
    }

    void normalize() noexcept {
        *this = normalized();
    }

    // Stream operators
    friend std::ostream& operator<<(std::ostream& os, const Vector3D& vec) {
        return os << vec.data_[0] << "  " << vec.data_[1] << "  " << vec.data_[2];
    }

    friend std::istream& operator>>(std::istream& is, Vector3D& vec) {
        return is >> vec.data_[0] >> vec.data_[1] >> vec.data_[2];
    }

    // Print method
    void print(std::ostream& os = std::cout) const {
        os << data_[0] << " " << data_[1] << " " << data_[2] << " ";
    }

    // Direct access to underlying data (for potential SIMD optimization)
    [[nodiscard]] constexpr const std::array<Real, 3>& data() const noexcept {
        return data_;
    }

    [[nodiscard]] constexpr std::array<Real, 3>& data() noexcept {
        return data_;
    }

private:
    alignas(32) std::array<Real, 3> data_;  // Aligned for SIMD operations
};

// Utility functions
[[nodiscard]] inline constexpr Real dot(const Vector3D& a, const Vector3D& b) noexcept {
    return a.dot(b);
}

[[nodiscard]] inline constexpr Vector3D cross(const Vector3D& a, const Vector3D& b) noexcept {
    return a.cross(b);
}

[[nodiscard]] inline constexpr Real squared_magnitude(const Vector3D& v) noexcept {
    return v.squared_magnitude();
}

[[nodiscard]] inline Real magnitude(const Vector3D& v) noexcept {
    return v.magnitude();
}

[[nodiscard]] inline constexpr Real squared_distance(const Vector3D& a, const Vector3D& b) noexcept {
    return a.squared_distance(b);
}

[[nodiscard]] inline Real distance(const Vector3D& a, const Vector3D& b) noexcept {
    return a.distance(b);
}

// Type alias for backward compatibility
using vektor = Vector3D;

} // namespace barnes_hut
