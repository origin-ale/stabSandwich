set term png size 850, 600 font ",16"
N = "100"
M = "12"

set output "output/deep_switch_optimization_graph_total.png"

set title "CAMPS-PP computation times, N×(N+12) circuit" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

set xrange [*:12.5]
set key bottom

plot 'output/deep_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "N=20, switch at t=32-c",\
'output/deep_switch_optimization_avgs_50.txt' using 1:2:3 with yerrorlines title "N=50, switch at t=62-c",\
'output/deep_switch_optimization_avgs_100.txt' using 1:2:3 with yerrorlines title "N=100, switch at t=112-c"