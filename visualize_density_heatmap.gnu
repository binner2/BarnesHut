#!/usr/bin/env gnuplot
#
# Modern Gnuplot 6 - Density Heatmap with 2D Histogram
# Shows particle concentration regions
#

set terminal pngcairo enhanced font "Arial,12" size 1600,1200 background rgb 'black'
set output 'particles_density_heatmap.png'

# Gnuplot 6 hexagonal binning for density estimation
set view map
set size ratio -1

# Modern inferno palette for heatmap
set palette defined (0 '#000004', 0.2 '#420a68', 0.4 '#932667', \
                     0.6 '#dd513a', 0.8 '#fca50a', 1 '#fcffa4')

# Title and labels
set title "Particle Density Heatmap (XY Projection)" \
    textcolor rgb 'white' font "Arial,16"
set xlabel "X Position" textcolor rgb 'white' font "Arial,14"
set ylabel "Y Position" textcolor rgb 'white' font "Arial,14"

# Border and tics
set border linecolor rgb 'white' linewidth 2
set xtics textcolor rgb 'white'
set ytics textcolor rgb 'white'

# Color box settings
set colorbox
set cblabel "Particle Count" textcolor rgb 'white'
set cbtics textcolor rgb 'white'

# Grid settings
set grid linecolor rgb '#404040'

# Use dgrid3d for density interpolation
set dgrid3d 50,50 qnorm 2
set pm3d map interpolate 0,0

# Plot with smoothed density
splot 'test.dat' skip 4 using 2:3:(1) with pm3d notitle

# Alternative: Point density visualization
# unset dgrid3d
# set style fill transparent solid 0.3
# plot 'test.dat' skip 4 using 2:3 with points pointtype 7 pointsize 1 \
#     linecolor rgb '#fde724' notitle
