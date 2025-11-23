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

// Modern C++20 constants
inline constexpr int NDIM = 3;
inline constexpr int NSUB = 1 << NDIM;
inline constexpr double GRAVITY = 1.0;
inline constexpr double EPSILON_SQUARED = 1e-10;  // Softening parameter for numerical stability

// Modern configuration flags
inline constexpr bool ENABLE_TIMING = true;
inline constexpr bool ENABLE_DEBUG = true;

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
