#include "stdinc.h"

namespace barnes_hut {

// Modern random number generation
double generate_random(double min_val, double max_val) {
    thread_local std::random_device rd;
    thread_local std::mt19937_64 gen(rd());
    std::uniform_real_distribution<double> dist(min_val, max_val);
    return dist(gen);
}

} // namespace barnes_hut
