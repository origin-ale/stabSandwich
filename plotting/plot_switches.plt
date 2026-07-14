set term png size 850, 600 font ",16"
N = "100"
M = "12"

set output "output/varM_switch_optimization_graph_total.png"

set title "CAMPS-PP computation times, 20×(20+M) circuit" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

set xrange [*:14.5]
set yrange [*:30]
set key bottom

plot 'output/switch/M6_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "M=6",\
'output/switch/M9_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "M=9",\
'output/switch/deep_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "M=12",\
'output/switch/M15_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "M=15",\
'output/switch/M18_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "M=18"
# 'output/switch/deep_switch_optimization_avgs_50.txt' using 1:2:3 with yerrorlines title "N=50, switch at t=62-c",\
# 'output/switch/deep_switch_optimization_avgs_100.txt' using 1:2:3 with yerrorlines title "N=100, switch at t=112-c"