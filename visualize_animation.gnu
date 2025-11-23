#!/usr/bin/env gnuplot
#
# Modern Gnuplot 6 - Animation Script (GIF)
# Creates animated visualization of particle evolution
#

# Animated GIF output
set terminal gif animate delay 10 size 1200,1200 \
    background rgb 'black' optimize

set output 'particles_animation.gif'

# Styling
set border linecolor rgb 'white' linewidth 2
set grid linecolor rgb '#404040' linewidth 1

# Fixed ranges for consistent animation
stats 'test.dat' skip 4 using 2 nooutput
xmin = STATS_min
xmax = STATS_max
stats 'test.dat' skip 4 using 3 nooutput
ymin = STATS_min
ymax = STATS_max

set xrange [xmin:xmax]
set yrange [ymin:ymax]

# Labels
set xlabel "X Position" textcolor rgb 'white' font "Arial,14"
set ylabel "Y Position" textcolor rgb 'white' font "Arial,14"
set xtics textcolor rgb 'white'
set ytics textcolor rgb 'white'

# Color palette
set palette defined (0 '#440154', 0.5 '#35b779', 1 '#fde724')

# Size ratio
set size ratio -1

# Animation loop - assuming we have multiple timestep files
# This is a template - actual implementation would need timestep data
do for [t=0:100:5] {
    set title sprintf("Barnes-Hut Simulation - Time Step %d", t) \
        textcolor rgb 'white' font "Arial,16"

    # In a real scenario, we'd load different files or sections
    # For now, we rotate the view for demonstration
    angle = t * 3.6
    set view 60, angle

    plot 'test.dat' skip 4 using 2:3:1 with points \
        pointtype 7 pointsize 1.0 palette notitle
}

# Note: For real time-evolution animation, you would need to:
# 1. Output particle positions at each timestep to separate files
# 2. Loop through these files in the animation
# Example: do for [i=0:num_steps] { plot sprintf('step_%04d.dat', i) ... }
