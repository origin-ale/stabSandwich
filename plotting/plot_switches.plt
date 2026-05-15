set term png size 850, 600 font ",16"
N = "60"

set output "output/switch_optimization_graph_".N.".png"

set title "Computation times, ".N."-qubit square circuit" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

set xrange [0:12.5]

plot 'output/switch_optimization_avgs_'.N.'.txt' using 1:2:3 with yerrorlines title "CAMPS-PP, switch at t=".N."-c"