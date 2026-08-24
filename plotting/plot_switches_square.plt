set term png size 850, 600 font ",16"

set output "output/varM_switch_optimization_graph_total.png"

set title "CAMPS-PP computation times, NxN circuit, switch at N-c" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

set xrange [*:14.5]
set yrange [*:30]
set key bottom

plot 'output/carlos_switch_square/M0_switch_optimization_avgs_20.txt' using 1:2:3 with yerrorlines title "N=20",\
'output/carlos_switch_square/M0_switch_optimization_avgs_50.txt' using 1:2:3 with yerrorlines title "N=50",\
'output/carlos_switch_square/M0_switch_optimization_avgs_80.txt' using 1:2:3 with yerrorlines title "N=80",\
'output/carlos_switch_square/M0_switch_optimization_avgs_100.txt' using 1:2:3 with yerrorlines title "N=100",\
