#!/usr/bin/env gnuplot
#
# Modern Gnuplot 6 - Velocity Field Visualization
# Vector field showing particle velocities
#

set terminal pngcairo enhanced font "Arial,12" size 1600,1200 background rgb 'black'
set output 'particles_velocity_field.png'

# Styling
set border linecolor rgb 'white' linewidth 2
set grid linecolor rgb '#404040' linewidth 1

# Title and labels
set title "Velocity Field (XY Projection)" \
    textcolor rgb 'white' font "Arial,16"
set xlabel "X Position" textcolor rgb 'white' font "Arial,14"
set ylabel "Y Position" textcolor rgb 'white' font "Arial,14"

# Tics
set xtics textcolor rgb 'white'
set ytics textcolor rgb 'white'

# Color palette based on velocity magnitude
set palette defined (0 '#1f77b4', 0.5 '#ff7f0e', 1 '#d62728')

# Key/Legend
set key textcolor rgb 'white'

# Size ratio
set size ratio -1

# Arrow style for vectors
set style arrow 1 head filled size screen 0.01,15 linewidth 1.5

# Calculate velocity magnitude and plot vectors
# Data format: mass x y z vx vy vz
# Sample every Nth particle to avoid clutter
plot 'test.dat' skip 4 every 10 using 2:3:(0.05*$5):(0.05*$6):(sqrt($5**2+$6**2)) \
    with vectors arrowstyle 1 palette title "Velocity Vectors", \
    '' skip 4 every 10 using 2:3 with points \
    pointtype 7 pointsize 0.5 linecolor rgb 'white' title "Particles"
