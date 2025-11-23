#!/usr/bin/env gnuplot
#
# Modern Gnuplot 6 - 2D Projections with Density Visualization
# XY, XZ, YZ projections with hexagonal binning
#

set terminal pngcairo enhanced font "Arial,11" size 2400,800 background rgb 'black'
set output 'particles_2d_projections.png'

# Multi-plot layout
set multiplot layout 1,3 title "Barnes-Hut Simulation - 2D Projections" \
    textcolor rgb 'white' font "Arial,18"

# Modern viridis-like palette
set palette defined (0 '#440154', 0.25 '#31688e', 0.5 '#35b779', 0.75 '#fde724', 1 '#ffffff')

# Common settings
set grid linecolor rgb '#404040' linewidth 1
set border linecolor rgb 'white' linewidth 2
set style fill transparent solid 0.5

#===============================================
# XY Projection
#===============================================
set title "XY Projection" textcolor rgb 'white' font "Arial,14"
set xlabel "X Position" textcolor rgb 'white'
set ylabel "Y Position" textcolor rgb 'white'
set xtics textcolor rgb 'white'
set ytics textcolor rgb 'white'
set size square

# Density plot with points colored by mass
plot 'test.dat' skip 4 using 2:3:1 with points \
    pointtype 7 pointsize 0.8 palette notitle

#===============================================
# XZ Projection
#===============================================
set title "XZ Projection" textcolor rgb 'white' font "Arial,14"
set xlabel "X Position" textcolor rgb 'white'
set ylabel "Z Position" textcolor rgb 'white'

plot 'test.dat' skip 4 using 2:4:1 with points \
    pointtype 7 pointsize 0.8 palette notitle

#===============================================
# YZ Projection
#===============================================
set title "YZ Projection" textcolor rgb 'white' font "Arial,14"
set xlabel "Y Position" textcolor rgb 'white'
set ylabel "Z Position" textcolor rgb 'white'

plot 'test.dat' skip 4 using 3:4:1 with points \
    pointtype 7 pointsize 0.8 palette notitle

unset multiplot
