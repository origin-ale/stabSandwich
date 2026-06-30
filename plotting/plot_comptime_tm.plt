set term png size 850, 600 font ",16"
set output "output/CompTimes_tm_graph.png"

set linetype 3 lc rgb "dark-red"

set title "N-qubit N/2-layer XXZ circuit" font ",24"
set logscale xy
set format y "%1g"
set xlabel "N"
set ylabel "Time (s)"

set xrange [1.9:20]

set key bottom

input = 'output/comptime_tm_cliff.txt'

plot input using 1:2 with linespoints title "CAMPS-PP",\
input using 1:3 with linespoints title "CAMPS",\
input using 1:4 with linespoints title "MPS",\
0.01*x**3 with lines lc rgb "gray" title "N^3"
