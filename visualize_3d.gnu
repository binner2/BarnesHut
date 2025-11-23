#!/usr/bin/env gnuplot
#
# Modern Gnuplot 6 Visualization for Barnes-Hut N-Body Simulation
# 3D Scatter Plot with Advanced Features
#

# Gnuplot 6 modern settings
set terminal pngcairo enhanced font "Arial,12" size 1920,1080 background rgb 'black'
set output 'particles_3d.png'

# Modern color palette with transparency
set palette defined (0 '#440154', 1 '#31688e', 2 '#35b779', 3 '#fde724')
set palette maxcolors 256

# Enhanced 3D settings
set view 60, 30, 1, 1
set xyplane at 0
set grid xtics ytics ztics linewidth 1.5 linecolor rgb '#404040'

# Axis settings with modern styling
set xlabel "X Position" textcolor rgb 'white' font "Arial,14"
set ylabel "Y Position" textcolor rgb 'white' font "Arial,14"
set zlabel "Z Position" textcolor rgb 'white' font "Arial,14"
set border linecolor rgb 'white' linewidth 2

# Title with metadata
set title "Barnes-Hut N-Body Simulation - 3D Particle Distribution" \
    textcolor rgb 'white' font "Arial,16"

# Tic marks styling
set xtics textcolor rgb 'white'
set ytics textcolor rgb 'white'
set ztics textcolor rgb 'white'

# Key/Legend
set key outside right textcolor rgb 'white'

# Point styling - modern approach with variable size
# Using pointtype 7 (filled circle) with size based on mass
set style fill transparent solid 0.6

# Data file format: mass x y z vx vy vz
# Skip first 4 lines (config data)
set datafile commentschars "#"

# 3D scatter plot with color mapping based on z-coordinate
splot 'test.dat' skip 4 using 2:3:4:4 with points \
    pointtype 7 pointsize 0.5 \
    palette notitle, \
    '' skip 4 using 2:3:4:(sprintf("%.0f", $1)) every 100 with labels \
    textcolor rgb 'yellow' font "Arial,8" offset 1,1 notitle
