set term png size 850, 600 font ",16"
set output "output/CompTimes_prob_tm_graph.png"

set linetype 3 lc rgb "dark-red"

set title "46-qubit 23-layer XXZ circuit" font ",24"
set logscale y
set format y "%1g"
set xlabel "magic probability"
set ylabel "Time (s)"

set xrange [0:*]

set key bottom right

input = 'output/comptime_prob_tm.txt'

plot input using 1:2:3 with yerrorlines title "CAMPS-PP",\
input using 1:4:5 with yerrorlines title "CAMPS",\
input using 1:6:7 with yerrorlines title "PP"
# input using 1:6:7 with yerrorlines title "MPS (thl 10^{-9})",\
