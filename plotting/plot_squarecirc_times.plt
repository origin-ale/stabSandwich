set term png size 850, 600 font ",16"
set output "output/CompTimes_squarecirc_graph.png"

set title "Computation times, N-qubit N-rotation circuit" font ",24"
set logscale xy
set format y "%.1e"
set xlabel "N"
set ylabel "Time (s)"

set key bottom

plot 'output/carlos_squarecirc3_v2/comptimes_3_squarecirc_avgs.txt' using 1:2:3 with yerrorlines title "Exact CAMPS",\
'output/carlos_squarecirc3_v2/comptimes_3_squarecirc_avgs.txt' using 1:4:5 with yerrorlines title "Exact PP",\
'output/carlos_squarecirc3_v2/comptimes_3_squarecirc_avgs.txt' using 1:6:7 with yerrorlines title "Exact MPS" lc rgb "red"