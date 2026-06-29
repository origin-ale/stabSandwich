set term png size 850, 600 font ",16"
set output "output/CompTimes_tm_graph.png"

set linetype 3 lc rgb "dark-red"

set title "N-qubit N/2-layer XXZ circuit" font ",24"
set logscale xy
set format y "%.1e"
set xlabel "N"
set ylabel "Time (s)"

set xrange [1.9:11]

set key bottom

plot 'output/comptime_tm.txt' using 1:2 with linespoints title "CAMPS-PP",\
'output/comptime_tm.txt' using 1:3 with linespoints title "CAMPS",\
'output/comptime_tm.txt' using 1:4 with linespoints title "MPS",\
0.001*x**4 with lines lc rgb "gray" title "x^4"
