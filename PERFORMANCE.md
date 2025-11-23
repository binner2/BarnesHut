# Performance Optimization Guide

## 🎯 Current Optimizations (Already Implemented)

### 1. CPU Parallelization (OpenMP)
**Status**: ✅ Implemented

```cpp
// Force calculation parallelized
#pragma omp parallel for schedule(dynamic)
for (Index i = 0; i < particles_.size(); ++i) {
    // Calculate forces for particle i
}

// Integration parallelized
#pragma omp parallel for
for (Index i = 0; i < particles_.size(); ++i) {
    particles_[i].integrate(dt_);
}
```

**Speedup**: 3.5-6.5× on 8-core CPU

### 2. Memory Optimization
**Status**: ✅ Implemented

- **Cache-line alignment**: 64-byte aligned Node structures
- **Pre-allocation**: Node pool reserves 3N nodes upfront
- **Memory reuse**: Nodes recycled between timesteps
- **Smart pointers**: Zero overhead with compiler optimization

### 3. Numerical Stability
**Status**: ✅ Implemented

```cpp
constexpr double EPSILON_SQUARED = 1e-10;  // Softening parameter
const Real r_cubed = (r_squared + EPSILON_SQUARED) * std::sqrt(r_squared + EPSILON_SQUARED);
```

**Benefit**: Prevents division by zero, improves accuracy

### 4. Compiler Optimizations
**Status**: ✅ Implemented in CMakeLists.txt

```cmake
-O3                    # Maximum optimization
-march=native          # CPU-specific instructions (AVX2, FMA)
-mtune=native          # CPU-specific tuning
-ffast-math            # Aggressive FP optimizations
-funroll-loops         # Loop unrolling
```

---

## 🚀 Future Performance Improvements

### Phase 1: Advanced CPU Optimizations

#### 1.1. SIMD Vectorization (AVX2/AVX-512)
**Potential Speedup**: 2-4× additional

**Current**:
```cpp
// Scalar operations
force += -GRAVITY * m1 * m2 / r_cubed * r_vec;
```

**Optimized**:
```cpp
#include <immintrin.h>

// Process 4 particles at once with AVX2
__m256d x = _mm256_load_pd(&positions_x[i]);
__m256d y = _mm256_load_pd(&positions_y[i]);
__m256d z = _mm256_load_pd(&positions_z[i]);
// ... SIMD force calculation
```

**Implementation Effort**: Medium
**Files to modify**: `tree.cpp`, `vektor.h`

#### 1.2. Structure of Arrays (SoA) Layout
**Potential Speedup**: 1.5-2× with SIMD

**Current (AoS)**:
```cpp
struct Particle {
    Real mass_;
    Vector3D position_;  // [x, y, z]
    Vector3D velocity_;
    Vector3D force_;
};
std::vector<Particle> particles;  // Non-contiguous x, y, z
```

**Optimized (SoA)**:
```cpp
struct ParticleArray {
    std::vector<Real> mass;
    std::vector<Real> pos_x, pos_y, pos_z;  // Contiguous!
    std::vector<Real> vel_x, vel_y, vel_z;
    std::vector<Real> force_x, force_y, force_z;
};
```

**Benefits**:
- Better cache utilization
- SIMD-friendly memory layout
- Reduced cache misses

**Implementation Effort**: High (requires refactoring)
**Files to modify**: `particle.h`, `tree.cpp`, all force calculations

#### 1.3. Task-Based Parallelism (Intel TBB)
**Potential Speedup**: 1.2-1.5× over OpenMP

Replace OpenMP with Intel Threading Building Blocks for better load balancing:

```cpp
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

tbb::parallel_for(tbb::blocked_range<size_t>(0, N),
    [&](const tbb::blocked_range<size_t>& r) {
        for (size_t i = r.begin(); i != r.end(); ++i) {
            calculate_forces(i);
        }
    }
);
```

**Benefits**:
- Better work stealing
- Recursive parallelism
- Lower overhead for small tasks

**Implementation Effort**: Medium
**Dependencies**: Intel TBB library

#### 1.4. Prefetching
**Potential Speedup**: 1.1-1.3×

```cpp
void interact(Particle& particle, const Node* node) {
    // Prefetch child nodes for next iteration
    for (const auto& child : node->children) {
        if (child) {
            __builtin_prefetch(child.get(), 0, 3);
        }
    }

    // Actual computation
    if (is_well_separated(particle, *node)) {
        particle_cell_interaction(particle, *node);
    }
}
```

**Implementation Effort**: Low
**Files to modify**: `tree.cpp`

---

### Phase 2: GPU Acceleration (CUDA/HIP)

#### 2.1. GPU Force Calculation
**Potential Speedup**: 10-50× over single CPU core

**Architecture**:
```
CPU (Host)              GPU (Device)
│                       │
├─ Build tree          │
├─ Compute CMS         │
│                       │
├─ Copy tree to GPU ──→│
│                       ├─ Parallel force calc
│                       │  (1 thread per particle)
│                       │
│←── Copy forces back ─┤
│                       │
├─ Integrate           │
└─ Update positions    │
```

**CUDA Kernel Example**:
```cuda
__global__ void calculate_forces_kernel(
    const Particle* particles,
    const Node* nodes,
    Vector3D* forces,
    int N
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;

    Vector3D force{0.0};

    // Tree traversal on GPU
    interact_gpu(particles[idx], &nodes[0], force);

    forces[idx] = force;
}
```

**Implementation Effort**: High
**Files to create**: `tree_gpu.cu`, `particle_gpu.cuh`
**Dependencies**: CUDA Toolkit or ROCm/HIP

#### 2.2. GPU-Accelerated Tree Build
**Potential Speedup**: 5-10×

Current tree building is serial. GPU version can parallelize:

```cuda
// Parallel Morton code generation
__global__ void compute_morton_codes(
    const Vector3D* positions,
    uint64_t* morton_codes,
    int N
);

// Parallel radix sort
thrust::sort_by_key(morton_codes, morton_codes + N, indices);

// Parallel tree construction
__global__ void build_tree_kernel(...);
```

**Implementation Effort**: Very High
**Algorithm**: Karras (2012) parallel BVH construction

#### 2.3. Multi-GPU Support
**Potential Speedup**: Near-linear scaling with GPU count

**Approach**:
- Domain decomposition
- Each GPU handles a spatial partition
- Halo exchange for boundary particles

**Implementation Effort**: Very High
**Dependencies**: NCCL (NVIDIA) or RCCL (AMD)

---

### Phase 3: Advanced Algorithmic Optimizations

#### 3.1. Fast Multipole Method (FMM)
**Potential Speedup**: O(N) vs O(N log N)

Barnes-Hut is O(N log N). FMM achieves O(N) by:
- Using multipole expansions up to higher orders
- Translation operators between cells
- Both upward and downward passes

**Implementation Effort**: Very High (research-level)
**Reference**: Greengard & Rokhlin (1987)

#### 3.2. Adaptive Time Stepping
**Potential Speedup**: 2-5× for realistic systems

Instead of fixed dt:
```cpp
// Individual timesteps per particle
Real dt_i = C * sqrt(eps / |a_i|);  // Courant condition

// Block time stepping (power-of-2 hierarchy)
if (time % dt_i == 0) {
    integrate(particle_i, dt_i);
}
```

**Benefits**:
- Tightly-bound systems use small dt
- Distant particles use large dt
- Overall fewer force calculations

**Implementation Effort**: High
**Files to modify**: `tree.cpp`, main loop

#### 3.3. Hierarchical Time Integration
**Potential Speedup**: 3-10× for multi-scale systems

```
Level 0: dt_0 = dt                    (fast-moving particles)
Level 1: dt_1 = 2×dt                  (medium)
Level 2: dt_2 = 4×dt                  (slow)
Level 3: dt_3 = 8×dt                  (very slow)
```

**Implementation Effort**: Very High
**Reference**: Springel (2005) - GADGET-2 code

#### 3.4. Tree Code + Particle Mesh (PM)
**Potential Speedup**: 2-3× for large-scale simulations

Hybrid approach:
- **Short-range**: Tree code (accurate)
- **Long-range**: FFT on grid (fast)

Split force:
```
F_total = F_tree(<r_split) + F_PM(>r_split)
```

**Implementation Effort**: Very High
**Dependencies**: FFTW library

---

## 📊 Performance Profiling

### Tools to Use

#### 1. perf (Linux)
```bash
perf record -g ./barnes_hut_sim data.dat 0.5 10
perf report
```

**Look for**:
- Hot functions (force calculation should dominate)
- Cache misses
- Branch mispredictions

#### 2. Intel VTune
```bash
vtune -collect hotspots ./barnes_hut_sim data.dat 0.5 10
```

**Metrics**:
- CPI (Cycles Per Instruction) - should be < 1.0
- L1/L2/L3 cache miss rates
- SIMD efficiency

#### 3. NVIDIA Nsight Compute (GPU)
```bash
ncu --set full ./barnes_hut_gpu data.dat 0.5 10
```

**Metrics**:
- SM efficiency (should be > 80%)
- Memory throughput
- Warp divergence

---

## 🎯 Optimization Priority Matrix

| Optimization | Speedup | Effort | Priority |
|-------------|---------|--------|----------|
| **Phase 1** | | | |
| SIMD (AVX2) | 2-4× | Medium | 🔥 High |
| SoA Layout | 1.5-2× | High | 🔥 High |
| Prefetching | 1.1-1.3× | Low | ⭐ Quick win |
| Intel TBB | 1.2-1.5× | Medium | ⭐ Quick win |
| **Phase 2** | | | |
| GPU Force Calc | 10-50× | High | 🚀 Huge |
| GPU Tree Build | 5-10× | Very High | 🚀 Huge |
| Multi-GPU | 2-4× | Very High | 🔬 Research |
| **Phase 3** | | | |
| FMM | 2-5× | Very High | 🔬 Research |
| Adaptive dt | 2-5× | High | 🔥 High |
| Hierarchical dt | 3-10× | Very High | 🔬 Research |
| Tree+PM | 2-3× | Very High | 🔬 Research |

---

## 💡 Quick Wins (Immediate Implementation)

### 1. Tune OpenMP Settings
```bash
export OMP_NUM_THREADS=8
export OMP_PROC_BIND=close
export OMP_PLACES=cores
```

### 2. Profile-Guided Optimization (PGO)
```bash
# Step 1: Build with profiling
g++ -fprofile-generate ...

# Step 2: Run with typical workload
./barnes_hut_sim typical_data.dat 0.5 10

# Step 3: Rebuild with profile data
g++ -fprofile-use ...
```

**Speedup**: 1.1-1.2× for free!

### 3. Link-Time Optimization (LTO)
```cmake
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
```

**Speedup**: 1.05-1.15×

### 4. Use jemalloc
```bash
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2 ./barnes_hut_sim ...
```

**Speedup**: 1.05-1.1× (better memory allocator)

---

## 🔬 Research-Level Optimizations

### 1. Machine Learning Acceleration
Use ML to predict whether to use multipole approximation:
```python
# Train neural network on (r, theta, accuracy) data
model.predict(r, theta) → use_approximation: bool
```

**Potential**: 1.2-1.5× by avoiding unnecessary tree traversals

### 2. Quantum Annealing (D-Wave)
For very large N, encode force calculation as QUBO problem.

**Potential**: Unknown (experimental)

### 3. Neuromorphic Computing
Map particles to neurons, forces to synapses.

**Potential**: Unknown (cutting-edge research)

---

## 📈 Expected Performance Trajectory

```
Current (2025):
├─ Serial: 1.0× baseline
├─ OpenMP (8 cores): 6× faster
└─ Optimized compiler flags: 1.5× additional

With Phase 1 (2026):
├─ + SIMD: 2.5× additional
├─ + SoA: 1.5× additional
└─ Total: ~20× vs serial baseline

With Phase 2 (2027):
├─ + Single GPU: 15× additional
├─ + Multi-GPU (4): 50× additional
└─ Total: ~1000× vs serial baseline

With Phase 3 (2028+):
├─ + FMM: 2× additional
├─ + Adaptive dt: 3× additional
└─ Total: ~6000× vs serial baseline
```

---

## 🏆 Summary

**Current State**: Modern C++20, OpenMP parallel, well-optimized

**Next Steps**:
1. ✅ Profile with perf/VTune
2. 🔥 Add SIMD vectorization (AVX2)
3. 🔥 Implement SoA memory layout
4. 🚀 Port to GPU (CUDA/HIP)

**Long-term Goal**: 1000× faster than 2004 code on modern hardware!
