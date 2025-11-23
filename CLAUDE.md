# CLAUDE.md - Barnes-Hut N-Body Simulation

## Project Overview

This is a **Barnes-Hut N-body simulation** implementation in C++ for computing gravitational/Coulomb forces between particles using the Fast Multipole Method (FMM). The code implements the Barnes-Hut algorithm with an octree data structure and uses Leapfrog integration for time-stepping.

**Original Author**: Burak İNNER
**Development Period**: 2004-2005
**Language**: C++
**Algorithm**: Barnes-Hut tree-based force calculation with multipole approximation

---

## Codebase Structure

### Core Components

```
BarnesHut/
├── BHtreetest.cpp      # Main program entry point and simulation loop
├── tree.h/cpp          # Barnes-Hut tree implementation (octree)
├── particle.h/cpp      # Particle class and Node structure definitions
├── vektor.h            # 3D vector class with mathematical operators
├── file.h/cpp          # File I/O for reading/writing simulation data
├── stdinc.h/cpp        # Standard includes and global constants
├── generate_data.cpp   # Utility to generate random test data
├── ncb                 # Binary data file (particle data)
└── .gitignore/.gitattributes
```

### File Descriptions

#### **BHtreetest.cpp**
- Main entry point for the Barnes-Hut simulation
- Handles command-line arguments: `./bh <filename> <TETA> <particleLEAF>`
- Implements the main simulation loop:
  1. Load particles into tree
  2. Upward pass (calculate center of mass)
  3. Calculate forces using tree approximation
  4. Empty tree for next iteration
  5. Output results
- Uses Leapfrog integration for time-stepping
- Contains timing/profiling code (controlled by `showtime` macro)

**Location References**:
- Main loop: `BHtreetest.cpp:83-158`
- Parameter validation: `BHtreetest.cpp:38-52`

#### **tree.h/cpp**
- Implements the Barnes-Hut octree data structure
- Key methods:
  - `loadparticles()`: Inserts particles into tree structure
  - `upwardpass()`: Computes center of mass for all nodes
  - `calculateforces()`: Computes forces using multipole approximation
  - `emptytree()`: Clears tree for next iteration
  - `addparticle()`: Recursively adds particle to correct octant
  - `whichchild()`: Determines which octant a particle belongs to
  - `leaf2node()`: Converts a leaf to internal node when capacity exceeded
- Uses STL `vector<Node>` for dynamic node allocation
- Implements well-separation criterion for multipole acceptance

**Key Algorithms**:
- Tree construction: `tree.cpp:loadparticles()`
- Center of mass calculation: `tree.cpp:upwardpass()`
- Force computation: `tree.cpp:calculateforces()`

#### **particle.h/cpp**
- Defines `particle` class with properties:
  - `data`: Mass or charge
  - `pos`: Position vector
  - `vel`: Velocity vector
  - `force`: Force accumulator
  - `Id`: Particle identifier
- Defines `Node` struct for tree nodes with:
  - `Type`: EMPTY, NODE, or LEAF
  - `geocenter`: Geometric center of octant
  - `rsize`: Size of octant (side length)
  - `masscenter`: Center of mass
  - `mass`: Total mass in subtree
  - `plist`: List of particles (if LEAF)
  - `child[NSUB]`: Pointers to 8 children (if NODE)
  - `parent`: Pointer to parent node

**Location References**:
- Node structure: `particle.h:16-40`
- Particle class: `particle.h:45-93`

#### **vektor.h**
- Implements 3D vector class with full operator overloading
- Supported operations:
  - Arithmetic: `+`, `-`, `*` (scalar and dot product), `/`
  - Compound assignment: `+=`, `-=`, `*=`, `/=`
  - Comparison: `==`, `!=`
  - Unary negation: `-v`
  - Array subscript: `v[i]` (0-indexed)
- Utility functions:
  - `square(v)`: Returns squared magnitude
  - `dist(v1, v2)`: Returns squared distance between vectors
  - I/O operators: `<<`, `>>`

**Note**: Uses inlined functions for performance

#### **file.h/cpp**
- File I/O operations for simulation data
- Key functions:
  - `readfilePart()`: Reads particle data from file
  - `putsnapshotfileForce()`: Writes force data to file
  - `putsnapshotfilePos()`: Writes position data to file
  - `generatedata()`: Generates random test data

**Data File Format**:
```
N                  # Number of particles
t_start            # Simulation start time
t_end              # Simulation end time
dt                 # Time step
mass x y z vx vy vz  # For each particle
mass x y z vx vy vz
...
```

#### **stdinc.h/cpp**
- Standard includes and global constants
- Key definitions:
  - `NDIM = 3`: Spatial dimensions
  - `NSUB = 8`: Number of children in octree (2^NDIM)
  - `gravity = 1`: Gravitational/Coulomb constant
  - `showtime`: Macro to enable timing output
  - `xdebug`: Macro to enable debugging output

#### **generate_data.cpp**
- Standalone utility to generate random particle data
- Creates test files with random masses, positions, and velocities
- Interactive program that prompts for:
  - Filename
  - Number of particles
  - Simulation start/end time
  - Time step
- Generates random data:
  - Mass: 5000-15000
  - Position: 0-10 in each dimension
  - Velocity: 0-100 in each dimension

---

## Key Algorithms

### 1. Barnes-Hut Tree Construction
- **Octree Subdivision**: Space is recursively divided into 8 octants
- **Particle Insertion**: Each particle is placed in the appropriate octant
- **Leaf Capacity**: When a leaf exceeds `particleLEAF` particles, it's split into 8 children

**Algorithm Flow**:
```
For each particle:
  1. Start at root
  2. Determine which octant contains particle (0-7)
  3. If octant is empty, create LEAF with particle
  4. If octant is LEAF:
     - If capacity not exceeded, add particle to LEAF
     - Otherwise, convert LEAF to NODE and redistribute particles
  5. If octant is NODE, recurse to that child
```

### 2. Upward Pass (Center of Mass Calculation)
- **Purpose**: Compute total mass and center of mass for each node
- **Method**: Post-order tree traversal (children before parent)
- **Formula**: `COM = Σ(m_i * pos_i) / Σ(m_i)`

### 3. Force Calculation
- **Multipole Approximation**: Distant groups of particles treated as single mass
- **Well-Separation Criterion**: `s/d < θ` where:
  - `s` = size of cell
  - `d` = distance to cell
  - `θ` = TETA parameter (typically 0.5-1.0)
- **Direct Calculation**: If not well-separated, calculate forces directly

**Algorithm Flow**:
```
For each particle:
  Walk tree from root:
    If node is well-separated:
      - Treat as single mass at center of mass
      - Calculate force contribution
    Else if node is LEAF:
      - Calculate direct forces with all particles in LEAF
    Else:
      - Recurse to children
```

### 4. Leapfrog Integration
- **Time-stepping scheme** for updating positions and velocities
- Second-order accurate, symplectic integrator
- Updates velocities at half-steps, positions at full steps

---

## Global Parameters

### TETA (θ - Theta Parameter)
- **Type**: `double`
- **Purpose**: Controls multipole approximation accuracy
- **Range**: Typically 0.0 - 2.0
- **Effect**:
  - Smaller values (0.3-0.5): More accurate, slower (more direct calculations)
  - Larger values (1.0-2.0): Faster, less accurate (more approximations)
- **Typical Value**: 0.5-1.0
- **Set via**: Command-line argument 2

### particleLEAF
- **Type**: `int`
- **Purpose**: Maximum number of particles allowed in a LEAF node
- **Range**: Typically 1-20
- **Effect**:
  - Smaller values: Deeper tree, more nodes, potentially faster force calculation
  - Larger values: Shallower tree, fewer nodes, more direct calculations in leaves
- **Typical Value**: 1-10
- **Set via**: Command-line argument 3

### Other Constants (stdinc.h)
- `NDIM = 3`: Spatial dimensions (DO NOT CHANGE without major code revision)
- `NSUB = 8`: Number of octree children (2^NDIM)
- `gravity = 1`: Force constant (can be adjusted for different simulations)

---

## Development Workflow

### Building the Project

**No Makefile is provided.** Manual compilation required.

#### Compile Main Simulation:
```bash
g++ -o bh BHtreetest.cpp tree.cpp particle.cpp file.cpp stdinc.cpp -O2 -Wall
```

#### Compile Data Generator:
```bash
g++ -o generate generate_data.cpp -O2 -Wall
```

#### Recommended Flags:
- `-O2` or `-O3`: Enable optimizations (critical for performance)
- `-Wall`: Enable all warnings
- `-g`: Add debug symbols (for debugging)
- `-std=c++11`: Use C++11 standard (recommended)

### Running the Simulation

```bash
# Generate test data (optional)
./generate

# Run simulation
./bh <datafile> <TETA> <particleLEAF>

# Example:
./bh filenamedata1.dat 0.7 5
```

**Arguments**:
1. `<datafile>`: Path to particle data file (e.g., `filenamedata1.dat`)
2. `<TETA>`: Theta parameter for multipole approximation (e.g., 0.7)
3. `<particleLEAF>`: Max particles per leaf node (e.g., 5)

### Expected Output

When `showtime` is defined:
- Load time for each step
- Upward pass time
- Force calculation time
- Total step time
- Statistics: direct force calculations, particle-cell interactions, nodes used

---

## Code Conventions

### Naming Conventions
- **Classes**: lowercase (e.g., `tree`, `particle`, `vektor`)
- **Functions**: camelCase (e.g., `loadparticles()`, `calculateforces()`)
- **Member variables**: camelCase (e.g., `geocenter`, `masscenter`)
- **Constants/Macros**: UPPERCASE (e.g., `NDIM`, `NSUB`, `TETA`)
- **Global variables**: UPPERCASE (e.g., `TETA`, `particleLEAF`)

### Code Style
- **Indentation**: Tabs (inconsistent in original code)
- **Braces**: K&R style (opening brace on same line)
- **Comments**: C-style `/* */` and C++ style `//`
- **Language**: Comments contain Turkish characters (original author's language)

### Memory Management
- Uses **dynamic node allocation** with STL `vector<Node>`
- Custom memory pool via `mynew()` and `Nodelist`
- **Important**: Tree is emptied and rebuilt each time step
- Particles stored as static array, allocated once

### Performance Considerations
- **Timing macros**: Use `#ifdef showtime` for profiling
- **Debug macros**: Use `#ifdef xdebug` for debugging
- **Inlined functions**: Vector operations are inlined for speed
- **STL containers**: `vector` used for dynamic node storage

---

## Important Implementation Details

### Tree Structure
- **Octree**: Each internal node has exactly 8 children (2^3)
- **Node Types**:
  - `EMPTY`: Unoccupied octant
  - `NODE`: Internal node with children
  - `LEAF`: Leaf node containing particles
- **Dynamic allocation**: Nodes allocated from vector pool via `mynew()`
- **Tree lifecycle**: Constructed each timestep, then emptied

### Particle-Tree Interaction
- Particles maintain reference to parent node via `parent` pointer
- LEAF nodes maintain `vector<particle*>` of contained particles
- Tree traversal used for force calculation, not particle storage

### Coordinate System
- **3D Cartesian**: Standard (x, y, z) coordinates
- **Indexing**: 0-based for vector elements
- **Octant numbering**: 0-7 based on binary representation of (x, y, z) quadrant

### Force Calculation
- **Coulomb/Gravity**: Same form (1/r^2), controlled by `gravity` constant
- **Softening**: No softening parameter (collision handling may be needed)
- **Units**: Dimensionless (user must ensure consistent units)

---

## Common Tasks for AI Assistants

### Adding New Features

#### 1. Adding a Softening Parameter
- **Where**: Modify `tree::directforcalc()` in `tree.cpp`
- **How**: Add `epsilon^2` to distance calculation: `r^2 + epsilon^2`
- **Also update**: `particle.h` to add softening member variable

#### 2. Implementing Periodic Boundaries
- **Where**: Modify `tree::whichchild()` and distance calculations
- **How**: Add modulo operations for wrapping coordinates
- **Files**: `tree.cpp`, `vektor.h` (modify `dist()` function)

#### 3. Adding Energy Calculation
- **Where**: Add method to `tree` class
- **How**: Traverse tree, accumulating potential energy using multipole approximation
- **Files**: `tree.h`, `tree.cpp`

#### 4. Supporting Different Force Laws
- **Where**: Modify `tree::directforcalc()` and `tree::particlecell()`
- **How**: Make force calculation a virtual function or template
- **Files**: `tree.cpp`

### Debugging Common Issues

#### Issue: Segmentation Fault
**Likely causes**:
1. NULL pointer dereference in tree traversal
2. Out-of-bounds vector access
3. Accessing cleared nodes

**Check**:
- Enable `xdebug` macro
- Verify `currentNode < Nodelist.size()`
- Check `parent` and `child` pointers before dereferencing

#### Issue: Incorrect Forces
**Likely causes**:
1. Wrong TETA parameter
2. Tree not properly constructed
3. Center of mass calculation error

**Check**:
- Print tree structure with `displaytree()`
- Verify upward pass completes before force calculation
- Compare with direct N^2 calculation for small N

#### Issue: Poor Performance
**Likely causes**:
1. TETA too small (too many direct calculations)
2. particleLEAF too large (not subdividing enough)
3. Compiler optimizations disabled

**Solutions**:
- Increase TETA to 0.7-1.0
- Decrease particleLEAF to 1-5
- Compile with `-O2` or `-O3`

### Modernization Opportunities

If updating this code to modern C++:
1. **Use smart pointers**: Replace raw `Node*` with `std::unique_ptr<Node>`
2. **Use std::array**: Replace `double element[3]` with `std::array<double, 3>`
3. **Use constexpr**: Make `NDIM`, `NSUB` constexpr instead of macros
4. **Add const correctness**: Many methods should be `const`
5. **Use nullptr**: Replace `NULL` with `nullptr`
6. **Add namespaces**: Wrap code in namespace to avoid global pollution
7. **Use enum class**: Replace `enum _type` with `enum class NodeType`
8. **Add CMake**: Create `CMakeLists.txt` for modern build system
9. **Remove timing macros**: Use proper profiling tools or runtime flags
10. **Fix Turkish encoding**: Convert comments to UTF-8 or English

---

## Testing Strategy

### Unit Testing Recommendations
1. **Vector operations**: Test all operators in `vektor.h`
2. **Tree construction**: Verify octant assignment for known positions
3. **Center of mass**: Check COM calculations for simple configurations
4. **Force calculation**: Compare with direct O(N^2) for small systems

### Integration Testing
1. **Two-body problem**: Should give exact forces (set TETA=0 for no approximation)
2. **Symmetric configurations**: Forces should be symmetric
3. **Energy conservation**: Total energy should be conserved (with symplectic integrator)

### Performance Testing
1. **Scaling**: Test with N = 100, 1000, 10000 particles
2. **TETA sensitivity**: Vary TETA and measure time vs accuracy
3. **particleLEAF tuning**: Find optimal value for target N

---

## File Format Reference

### Input Data File Format
```
N                          # Line 1: Number of particles (integer)
t_start                    # Line 2: Start time (double)
t_end                      # Line 3: End time (double)
dt                         # Line 4: Time step (double)
mass1 x1 y1 z1 vx1 vy1 vz1 # Line 5+: Particle data (7 doubles per line)
mass2 x2 y2 z2 vx2 vy2 vz2
...
```

**Example**:
```
100
0.0
10.0
0.01
5000.0 5.3 2.1 8.9 10.0 -5.0 3.2
7500.0 1.2 9.8 3.4 -8.0 12.0 -1.5
...
```

### Output Files
- Created by `putsnapshotfileForce()` and `putsnapshotfilePos()`
- Contains force or position vectors for all particles
- Format determined by implementation in `file.cpp`
- Used for visualization and analysis

---

## Dependencies

### Standard Library Headers
- `<iostream>`: Console I/O
- `<fstream>`: File I/O
- `<cmath>`: Mathematical functions (sqrt, etc.)
- `<cstdlib>`: Standard utilities (exit, atof, atoi)
- `<ctime>`: Timing functions (clock, CLOCKS_PER_SEC)
- `<iomanip>`: I/O manipulators
- `<vector>`: STL vector container

### External Dependencies
**None** - This is a standalone C++ project with no external libraries.

---

## Performance Characteristics

### Computational Complexity
- **Tree construction**: O(N log N)
- **Force calculation**: O(N log N) with Barnes-Hut
- **Direct calculation**: O(N^2) (used as fallback)
- **Space complexity**: O(N) for particles, O(N) for tree nodes

### Typical Performance
For N = 10,000 particles:
- Tree construction: ~0.1-0.5 seconds
- Force calculation: ~0.5-2.0 seconds per step
- Total: ~1-3 seconds per timestep (hardware dependent)

Performance scales well with optimizations enabled (`-O2`, `-O3`).

---

## Known Issues and Limitations

1. **No Makefile**: Manual compilation required
2. **Global variables**: `TETA` and `particleLEAF` are global (not thread-safe)
3. **No error handling**: File I/O assumes files exist and are valid
4. **Turkish comments**: Some comments in original author's language
5. **Encoding issues**: Non-ASCII characters may not display correctly
6. **No softening**: Particles can experience infinite forces at r=0
7. **Single-threaded**: No parallelization (could benefit from OpenMP)
8. **Memory management**: Tree rebuilt each step (could be optimized)

---

## Quick Reference for AI Assistants

### When modifying this code:
1. **Always** enable compiler warnings (`-Wall`)
2. **Test** with small N before running large simulations
3. **Verify** energy conservation for physical accuracy
4. **Profile** before optimizing (use timing output)
5. **Preserve** the tree construction → upward pass → force calculation order
6. **Never** access tree nodes after `emptytree()` is called
7. **Be careful** with global variables `TETA` and `particleLEAF`
8. **Remember** that `NSUB = 8` is hardcoded for 3D (changing NDIM requires code changes)

### Common gotchas:
- Vector indices are 0-based: `v[0]`, `v[1]`, `v[2]`
- `dist()` returns **squared** distance (no sqrt)
- Node children indexed 0-7 (binary representation of octant)
- Force accumulation, not assignment (use `+=` not `=`)
- Tree is rebuilt each timestep (not persistent)

### Files to modify for:
- **Algorithm changes**: `tree.cpp`
- **Data structures**: `particle.h`, `tree.h`
- **I/O format**: `file.cpp`
- **Constants**: `stdinc.h`
- **Integration scheme**: `particle.cpp` (currently minimal)
- **Main loop**: `BHtreetest.cpp`

---

**Last Updated**: 2025-11-23
**For**: AI assistants working on this Barnes-Hut N-body simulation code
