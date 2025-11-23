#!/bin/bash
#
# Master visualization script for Barnes-Hut simulation
# Generates all plots using Gnuplot 6
#

echo "=== Barnes-Hut Visualization Suite ==="
echo ""

# Check if gnuplot is installed
if ! command -v gnuplot &> /dev/null; then
    echo "Error: gnuplot not found. Please install gnuplot 6.0+"
    exit 1
fi

# Check gnuplot version
GNUPLOT_VERSION=$(gnuplot --version | grep -oP '\d+\.\d+' | head -1)
echo "Using gnuplot version: $GNUPLOT_VERSION"
echo ""

# Check if data file exists
if [ ! -f "test.dat" ]; then
    echo "Generating test data..."
    ./generate_data test.dat 1000 0.0 1.0 0.01
    echo ""
fi

echo "Generating visualizations..."
echo ""

# 1. 3D Scatter Plot
echo "[1/5] Creating 3D scatter plot..."
gnuplot visualize_3d.gnu
if [ $? -eq 0 ]; then
    echo "  ✓ particles_3d.png created"
else
    echo "  ✗ Failed to create 3D plot"
fi

# 2. 2D Projections
echo "[2/5] Creating 2D projections..."
gnuplot visualize_2d_projections.gnu
if [ $? -eq 0 ]; then
    echo "  ✓ particles_2d_projections.png created"
else
    echo "  ✗ Failed to create 2D projections"
fi

# 3. Density Heatmap
echo "[3/5] Creating density heatmap..."
gnuplot visualize_density_heatmap.gnu
if [ $? -eq 0 ]; then
    echo "  ✓ particles_density_heatmap.png created"
else
    echo "  ✗ Failed to create density heatmap"
fi

# 4. Velocity Field
echo "[4/5] Creating velocity field..."
gnuplot visualize_velocity_field.gnu
if [ $? -eq 0 ]; then
    echo "  ✓ particles_velocity_field.png created"
else
    echo "  ✗ Failed to create velocity field"
fi

# 5. Animation (optional - can be slow)
read -p "[5/5] Create animation? (y/n): " -n 1 -r
echo ""
if [[ $REPLY =~ ^[Yy]$ ]]; then
    echo "Creating animation (this may take a while)..."
    gnuplot visualize_animation.gnu
    if [ $? -eq 0 ]; then
        echo "  ✓ particles_animation.gif created"
    else
        echo "  ✗ Failed to create animation"
    fi
else
    echo "  ⊗ Animation skipped"
fi

echo ""
echo "=== Visualization Complete ==="
echo "Generated files:"
ls -lh particles_*.png particles_*.gif 2>/dev/null | awk '{print "  " $9 " (" $5 ")"}'
