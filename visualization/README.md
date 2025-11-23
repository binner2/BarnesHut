# CUDA-Accelerated Visualization for Barnes-Hut N-Body Simulation

Modern real-time visualization system combining CUDA compute with OpenGL 4.5+ rendering.

## Features

- **CUDA-OpenGL Interop**: Zero-copy data transfer between CUDA and OpenGL
- **GPU-Accelerated Rendering**: Millions of particles rendered in real-time
- **Dynamic Color Mapping**: Visualize velocity or force magnitude with smooth color gradients
- **Interactive Camera**: Orbit, pan, and zoom with mouse controls
- **Modern Shaders**: Point sprites with soft edges, glow effects, and distance attenuation
- **Real-time Controls**: Pause, step, adjust rendering parameters on-the-fly

## Architecture

```
┌─────────────────┐
│  Particle Data  │
└────────┬────────┘
         │
         ├──> CPU: Barnes-Hut Tree Construction
         │         (OpenMP parallelized)
         │
         └──> GPU: CUDA Kernel
              ├─> Color Mapping
              ├─> Position Transform
              └─> OpenGL VBO (Direct Write)
                  │
                  └──> OpenGL Rendering
                       ├─> Vertex Shader
                       ├─> Fragment Shader
                       └─> Screen Output
```

## Dependencies

### Required

1. **CUDA Toolkit 11.0+**
   ```bash
   # Check if installed
   nvcc --version

   # Install on Ubuntu
   wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
   sudo dpkg -i cuda-keyring_1.1-1_all.deb
   sudo apt-get update
   sudo apt-get install cuda-toolkit-12-0
   ```

2. **OpenGL 4.5+** (usually provided by GPU drivers)

3. **GLFW 3**
   ```bash
   sudo apt-get install libglfw3-dev
   ```

4. **GLM** (OpenGL Mathematics)
   ```bash
   sudo apt-get install libglm-dev
   ```

5. **GLAD** (OpenGL Function Loader)
   ```bash
   # Download GLAD
   cd visualization
   mkdir -p glad
   cd glad

   # Option 1: Use the web service
   # Visit https://glad.dav1d.de/
   # Settings:
   #   - Language: C/C++
   #   - Specification: OpenGL
   #   - gl: Version 4.5 (or higher)
   #   - Profile: Core
   #   - Extensions: None (or add specific ones)
   #   - Options: Check "Generate a loader"
   # Download the generated files and extract here

   # Option 2: Use glad generator (Python)
   pip install glad
   glad --api="gl:core=4.5" --out-path=. --generator=c

   # You should have:
   #   glad/glad.h
   #   glad/glad.c
   #   glad/KHR/khrplatform.h
   ```

### Optional

- **OpenMP** (for CPU parallelization)
  ```bash
  sudo apt-get install libomp-dev
  ```

## Building

```bash
# Create build directory
mkdir build
cd build

# Configure with CMake
cmake .. -DCMAKE_BUILD_TYPE=Release

# Build
make -j$(nproc)

# The executable will be: barnes_hut_visual
```

### CMake Options

```bash
# Specify CUDA architecture (optional)
cmake .. -DCMAKE_CUDA_ARCHITECTURES="75;80;86"

# Disable visualization if dependencies are missing
cmake .. -DENABLE_VISUALIZATION=OFF

# Debug build
cmake .. -DCMAKE_BUILD_TYPE=Debug
```

## Usage

```bash
./barnes_hut_visual <data_file> <theta> <particles_per_leaf>
```

### Example

```bash
# Generate test data
./generate_data
# Enter: 10000 particles, 0.0 start time, 10.0 end time, 0.01 time step

# Run visualization
./barnes_hut_visual data.dat 0.7 8
```

## Controls

### Keyboard

| Key | Action |
|-----|--------|
| `Space` | Pause/Resume simulation |
| `S` | Single step (when paused) |
| `R` | Reset camera |
| `V` | Toggle velocity/force coloring |
| `G` | Toggle glow effect |
| `+/-` | Increase/Decrease particle size |
| `[/]` | Decrease/Increase simulation speed |
| `H` | Toggle help |
| `ESC` | Exit |

### Mouse

| Action | Control |
|--------|---------|
| Rotate | Left drag |
| Pan | Right drag |
| Zoom | Mouse wheel |

## Performance Tips

1. **GPU Selection**: Ensure you're using the NVIDIA GPU
   ```bash
   # Check available GPUs
   nvidia-smi

   # If using laptop with integrated + discrete GPU
   export __NV_PRIME_RENDER_OFFLOAD=1
   export __GLX_VENDOR_LIBRARY_NAME=nvidia
   ./barnes_hut_visual data.dat 0.7 8
   ```

2. **Particle Count**: Start with smaller datasets (1K-10K particles) to test
   - Small (< 10K): Real-time with high quality
   - Medium (10K-100K): Real-time with optimizations
   - Large (> 100K): May need to reduce rendering quality

3. **CUDA Architecture**: Compile for your specific GPU
   ```bash
   # Find your GPU architecture
   nvidia-smi --query-gpu=compute_cap --format=csv

   # Add to CMake
   cmake .. -DCMAKE_CUDA_ARCHITECTURES="86"  # For RTX 3090/3080
   ```

4. **VSync**: Disable for maximum frame rate
   - Edit `visualization_main.cpp`: Set `window_config.vsync = false;`

## Troubleshooting

### "CUDA error: no CUDA-capable device is detected"
- Ensure NVIDIA drivers are installed: `nvidia-smi`
- Check CUDA installation: `nvcc --version`

### "Failed to load shaders"
- Ensure you're running from the build directory
- Check that `shaders/` directory exists with `.vert` and `.frag` files
- CMake should copy shaders automatically

### "Failed to initialize GLAD"
- Ensure GLAD files are in `visualization/glad/`
- Verify OpenGL 4.5+ support: `glxinfo | grep "OpenGL version"`

### "GLFW3 not found"
```bash
sudo apt-get update
sudo apt-get install libglfw3-dev
```

### "GLM not found"
```bash
sudo apt-get install libglm-dev
```

## Color Schemes

### Velocity Coloring (Default)
- **Blue**: Slow particles
- **Cyan**: Moderate speed
- **Green**: Medium-fast
- **Yellow**: Fast
- **Red**: Very fast

### Force Coloring
- **White**: Low force
- **Yellow**: Moderate force
- **Orange**: High force
- **Red**: Very high force
- **Magenta**: Extreme force

## Technical Details

### CUDA Kernels

- `prepare_vertices_kernel`: Converts particle data to OpenGL vertex format with color mapping
- `compute_velocity_stats_kernel`: Parallel reduction for velocity statistics

### OpenGL Shaders

- **Vertex Shader** (`particle.vert`): Transforms particles to clip space with distance-based size attenuation
- **Fragment Shader** (`particle.frag`): Renders circular point sprites with soft edges and optional glow

### Memory Layout

```cpp
struct VertexData {
    float x, y, z;      // Position (12 bytes)
    float r, g, b, a;   // Color RGBA (16 bytes)
};  // Total: 28 bytes per particle (aligned to 32 bytes)
```

### Interop Flow

```
CPU Particles → CUDA Memory → CUDA Kernel → OpenGL VBO (mapped) → GPU Render
     |                                              ↑
     └──────────────────────────────────────────────┘
              (No CPU-GPU copy overhead)
```

## Future Enhancements

- [ ] ImGui interface for runtime parameter adjustment
- [ ] Multiple rendering modes (points, billboards, volumetric)
- [ ] Bloom post-processing effect
- [ ] Recording to video
- [ ] Multi-GPU support
- [ ] Compute shader alternative to CUDA
- [ ] Vulkan backend

## License

Same as parent project (Barnes-Hut simulation).

## Credits

- **CUDA-OpenGL Interop**: NVIDIA CUDA Samples
- **Modern OpenGL**: LearnOpenGL (learnopengl.com)
- **Camera Controller**: Inspired by orbit camera implementations
