set term png size 850, 600 font ",16"
set output "output/CompTimes_mps_graph.png"

set linetype 3 lc rgb "dark-red"

set title "MPS computation times, N-qubit N-rotation circuit" font ",24"
set logscale xy
set format y "%.1e"
set xlabel "N"
set ylabel "Time (s)"

set key bottom

plot 'output/comptimes_squarecirc_mps_ex_avgs.txt' using 1:2:3 with yerrorlines title "Exact",\
'output/comptimes_squarecirc_mps_20_avgs.txt' using 1:2:3 with yerrorlines title "SVD cutoff 10^{-20}",\
'output/comptimes_squarecirc_mps_12_avgs.txt' using 1:2:3 with yerrorlines title "SVD cutoff 10^{-12}",\
'output/comptimes_squarecirc_mps_6_avgs.txt' using 1:2:3 with yerrorlines title "SVD cutoff 10^{-6}"