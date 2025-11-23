#pragma once

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <iomanip>
#include <string>
#include <memory>
#include <vector>
#include <array>
#include <optional>
#include <string_view>
#include <span>
#include <concepts>
#include <execution>
#include <algorithm>
#include <random>
#include <ranges>

// Modern C++20 constants
inline constexpr int NDIM = 3;
inline constexpr int NSUB = 1 << NDIM;  // 2^3 = 8 octants
inline constexpr double GRAVITY = 1.0;
inline constexpr double EPSILON_SQUARED = 1e-10;  // Softening parameter for numerical stability

// Modern configuration flags
inline constexpr bool ENABLE_TIMING = true;
inline constexpr bool ENABLE_DEBUG = true;

// Compile-time mathematical constants
namespace math {
    inline constexpr double PI = 3.14159265358979323846;
    inline constexpr double TWO_PI = 2.0 * PI;
    inline constexpr double HALF_PI = 0.5 * PI;
    inline constexpr double E = 2.71828182845904523536;
    inline constexpr double SQRT_2 = 1.41421356237309504880;
    inline constexpr double SQRT_3 = 1.73205080756887729352;

    // Compile-time power function
    template<typename T>
    constexpr T pow(T base, unsigned int exp) noexcept {
        return (exp == 0) ? 1 : base * pow(base, exp - 1);
    }

    // Compile-time square
    template<typename T>
    constexpr T square(T x) noexcept {
        return x * x;
    }

    // Compile-time cube
    template<typename T>
    constexpr T cube(T x) noexcept {
        return x * x * x;
    }

    // Fast inverse square root approximation (Quake III algorithm)
    // Using std::bit_cast (C++20) for safe type punning
    inline float fast_inv_sqrt(float x) noexcept {
        const float x_half = 0.5f * x;
        int i = std::bit_cast<int>(x);
        i = 0x5f3759df - (i >> 1);
        x = std::bit_cast<float>(i);
        x = x * (1.5f - x_half * x * x);  // One Newton iteration
        return x;
    }
}  // namespace math

namespace barnes_hut {

// Modern type aliases
using Real = double;
using Index = std::size_t;
using TimePoint = std::chrono::high_resolution_clock::time_point;
using Duration = std::chrono::duration<double>;

// High-precision timer class
class Timer {
public:
    Timer() : start_(std::chrono::high_resolution_clock::now()) {}

    void reset() { start_ = std::chrono::high_resolution_clock::now(); }

    [[nodiscard]] double elapsed() const {
        auto end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start_).count();
    }

private:
    TimePoint start_;
};

} // namespace barnes_hut
