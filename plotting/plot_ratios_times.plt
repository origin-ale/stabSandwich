set term png size 850, 600 font ",16"
N = "40"

set output "output/CompTimes_ratios_".N."_graph.png"

set title "Computation times, N=".N." t-rotation circuit" font ",24"
set xrange [0.65:1.25]
set logscale y
set format y "%.1e"
set xlabel "t/N"
set ylabel "Time (s)"

set key bottom

plot 'output/comptimes_ratios_'.N.'_avgs.txt' using 1:2:3 with yerrorlines title "Exact CAMPS",\
'output/comptimes_ratios_'.N.'_avgs.txt' using 1:4:5 with yerrorlines title "Exact PP"