# Building Barnes-Hut with CUDA Visualization

Quick start guide for building the Barnes-Hut simulation with modern CUDA-accelerated visualization.

## Prerequisites Check

```bash
# Check CUDA installation
nvcc --version  # Should show CUDA 11.0 or higher

# Check NVIDIA driver
nvidia-smi      # Should show your GPU

# Check OpenGL version
glxinfo | grep "OpenGL version"  # Should show 4.5 or higher
```

## Step-by-Step Build Instructions

### 1. Install Dependencies

#### Ubuntu/Debian

```bash
# Update package list
sudo apt-get update

# Install build essentials
sudo apt-get install -y build-essential cmake git

# Install CUDA Toolkit (if not already installed)
# Visit: https://developer.nvidia.com/cuda-downloads
# Or use:
sudo apt-get install -y nvidia-cuda-toolkit

# Install OpenGL and windowing libraries
sudo apt-get install -y libglfw3-dev libglm-dev libomp-dev

# Install OpenGL utilities (for glxinfo)
sudo apt-get install -y mesa-utils
```

#### Fedora/RHEL

```bash
sudo dnf install -y gcc-c++ cmake git
sudo dnf install -y cuda cuda-toolkit
sudo dnf install -y glfw-devel glm-devel libomp-devel
sudo dnf install -y glx-utils
```

#### Arch Linux

```bash
sudo pacman -S base-devel cmake git
sudo pacman -S cuda glfw-x11 glm openmp
sudo pacman -S mesa-utils
```

### 2. Setup GLAD (OpenGL Function Loader)

GLAD is required for loading OpenGL functions. Choose one method:

#### Method A: Using the provided script (Recommended)

```bash
cd visualization
./download_glad.sh

# If you have Python pip:
pip install glad
./download_glad.sh
```

#### Method B: Manual download

1. Visit https://glad.dav1d.de/
2. Configure:
   - Language: **C/C++**
   - Specification: **OpenGL**
   - gl: **Version 4.5** (or higher)
   - Profile: **Core**
   - Options: ✓ **Generate a loader**
3. Click **Generate** and download
4. Extract to `visualization/glad/`:
   ```bash
   cd visualization
   mkdir -p glad
   # Extract downloaded files here
   # You should have: glad/glad.h, glad/glad.c, glad/KHR/khrplatform.h
   ```

#### Method C: Using Python glad

```bash
pip install glad
cd visualization
mkdir -p glad
cd glad
glad --api="gl:core=4.5" --out-path=. --generator=c

# Reorganize files
mv glad/glad.h .
mv glad/KHR .
mv src/glad.c .
rm -rf glad src
```

### 3. Build the Project

```bash
# From project root directory
mkdir -p build
cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build (use all CPU cores)
make -j$(nproc)

# You should now have three executables:
# - barnes_hut_sim (original CLI version)
# - barnes_hut_visual (new visualization version)
# - generate_data (data generator)
```

### 4. Generate Test Data

```bash
# From build directory
./generate_data

# When prompted, enter:
# Filename: test_10k.dat
# Number of particles: 10000
# Start time: 0.0
# End time: 10.0
# Time step: 0.01
```

### 5. Run Visualization

```bash
# Syntax: ./barnes_hut_visual <data_file> <theta> <particles_per_leaf>
./barnes_hut_visual test_10k.dat 0.7 8
```

## Expected Output

```
=== Barnes-Hut N-Body Visualization ===

Loading simulation data...
Loaded 10000 particles
  Simulation time: 0 -> 10
  Time step: 0.01
  Theta: 0.7
  Particles per leaf: 8

Initializing visualization system...
Window initialized: 1920x1080
OpenGL Renderer initialized for 10000 particles
  OpenGL Version: 4.5.0 NVIDIA
  GLSL Version: 4.50 NVIDIA
  Renderer: NVIDIA GeForce RTX 3080
CUDA Renderer initialized for 10000 particles
OpenGL VBO registered with CUDA
Visualization initialized successfully

=== Keyboard Controls ===
  Space     : Pause/Resume simulation
  ...

Starting visualization loop...
```

## Troubleshooting

### Issue: "CUDA not found"

```bash
# Check CUDA installation
which nvcc
# If not found, add to PATH:
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

# Add to ~/.bashrc to make permanent
```

### Issue: "GLFW3 not found"

```bash
# Install GLFW
sudo apt-get install libglfw3-dev

# If still not found, try:
sudo apt-get install libglfw3
sudo apt-get install pkg-config
```

### Issue: "GLM not found"

```bash
# Install GLM (header-only library)
sudo apt-get install libglm-dev
```

### Issue: "Failed to load shaders"

```bash
# Make sure you're in the build directory
cd build

# Check if shaders directory exists
ls shaders/
# Should show: particle.vert  particle.frag

# If not, copy manually:
cp -r ../visualization/shaders .
```

### Issue: "no CUDA-capable device"

```bash
# Check NVIDIA driver
nvidia-smi

# If not working, reinstall driver:
sudo apt-get install --reinstall nvidia-driver-XXX  # Replace XXX with your version
```

### Issue: Black screen or no rendering

```bash
# Check OpenGL version
glxinfo | grep "OpenGL version"
# Must be 4.5 or higher

# For laptops with integrated + discrete GPU:
export __NV_PRIME_RENDER_OFFLOAD=1
export __GLX_VENDOR_LIBRARY_NAME=nvidia
./barnes_hut_visual test_10k.dat 0.7 8
```

### Issue: Slow performance

1. **Disable VSync** for maximum FPS:
   - Edit `visualization/visualization_main.cpp`
   - Change `window_config.vsync = false;`
   - Rebuild

2. **Check GPU usage**:
   ```bash
   nvidia-smi -l 1  # Monitor GPU usage
   ```

3. **Reduce particle count** for testing:
   ```bash
   ./generate_data  # Create smaller dataset (1000 particles)
   ```

4. **Compile for your GPU architecture**:
   ```bash
   # Check your GPU compute capability
   nvidia-smi --query-gpu=compute_cap --format=csv

   # Rebuild with specific architecture (e.g., 8.6 for RTX 30 series)
   cmake .. -DCMAKE_CUDA_ARCHITECTURES="86"
   make -j$(nproc)
   ```

## Build Options

### Disable Visualization

If you want to build only the CPU version without visualization:

```bash
cmake .. -DENABLE_VISUALIZATION=OFF
make -j$(nproc)
```

### Debug Build

For debugging with CUDA-GDB:

```bash
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j$(nproc)

# Run with cuda-gdb
cuda-gdb --args ./barnes_hut_visual test_10k.dat 0.7 8
```

### Specific CUDA Architectures

```bash
# For RTX 30 series (Ampere - compute capability 8.6)
cmake .. -DCMAKE_CUDA_ARCHITECTURES="86"

# For RTX 20 series (Turing - compute capability 7.5)
cmake .. -DCMAKE_CUDA_ARCHITECTURES="75"

# For multiple architectures
cmake .. -DCMAKE_CUDA_ARCHITECTURES="75;80;86"
```

## Performance Benchmarks

Expected performance on various hardware:

| GPU | Particles | FPS |
|-----|-----------|-----|
| RTX 4090 | 100K | 60+ |
| RTX 3080 | 100K | 50+ |
| RTX 3060 | 50K | 60+ |
| RTX 2060 | 10K | 60+ |
| GTX 1660 | 5K | 60+ |

*Note: FPS depends on theta value, particles_per_leaf, and rendering settings*

## Next Steps

After successful build:

1. **Read the controls**: Press `H` in the visualization window
2. **Experiment with settings**:
   - Try different theta values (0.5 - 1.0)
   - Adjust particle sizes with `+/-`
   - Toggle color modes with `V`
3. **Create larger simulations**:
   ```bash
   ./generate_data  # Try 100K particles
   ./barnes_hut_visual large_data.dat 0.9 16
   ```

## Additional Resources

- **Visualization README**: `visualization/README.md`
- **Main Project README**: `README.md`
- **CUDA Programming Guide**: https://docs.nvidia.com/cuda/
- **OpenGL Tutorial**: https://learnopengl.com/

## Getting Help

If you encounter issues:

1. Check this troubleshooting guide
2. Verify all dependencies are installed
3. Check `visualization/README.md` for more details
4. Open an issue on GitHub with:
   - Your GPU model (`nvidia-smi`)
   - CUDA version (`nvcc --version`)
   - OpenGL version (`glxinfo | grep "OpenGL version"`)
   - Full error message

## Summary

```bash
# Quick build (Ubuntu):
sudo apt-get install -y libglfw3-dev libglm-dev nvidia-cuda-toolkit
cd visualization && ./download_glad.sh && cd ..
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j$(nproc)
./generate_data
./barnes_hut_visual test.dat 0.7 8
```

Enjoy the visualization! 🚀
