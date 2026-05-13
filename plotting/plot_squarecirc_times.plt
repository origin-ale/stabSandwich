set term png size 850, 600 font ",16"
set output "output/CompTimes_squarecirc_graph.png"

set title "Computation times, N-qubit N-rotation circuit" font ",24"
set logscale xy
set format y "%.1e"
set xlabel "N"
set ylabel "Time (s)"

plot 'output/comptimes_squarecirc_avgs.txt' using 1:2:3 with yerrorlines title "Exact CAMPS",\
'output/comptimes_squarecirc_avgs.txt' using 1:4:5 with yerrorlines title "Exact PP"