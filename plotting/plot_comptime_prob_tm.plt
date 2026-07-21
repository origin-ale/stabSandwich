set term png size 850, 600 font ",16"
set output "output/CompTimes_prob_tm_graph.png"

set linetype 3 lc rgb "dark-red"

set title "46-qubit 23-layer XXZ circuit" font ",24"
set logscale y
set format y "%1g"
set xlabel "Magic probability"
set ylabel "Time (s)"

set xrange [0:*]
set yrange [0:1000]

set key top left

plot 'output/carlos_varseed64/comptime_prob_tm.txt' using 1:2:3 with yerrorlines title "CAMPS-PP, χ=64, thl=10^{-10}",\
'output/carlos_varseed_pp/comptime_prob_tm.txt' using 1:2:3 with yerrorlines title "PP, thl 10^{-10}"
# plot 'output/carlos_varseed/comptime_prob_tm.txt' using 1:2:3 with yerrorlines title "CAMPS-PP, χ=128, thl=10^{-10}",\
# 'output/carlos_varseed32/comptime_prob_tm.txt' using 1:2:3 with yerrorlines title "CAMPS-PP, χ=32, thl=10^{-10}",\
# input using 1:4:5 with yerrorlines title "CAMPS - search"\,
# input using 1:6:7 with yerrorlines title "MPS (thl 10^{-9})",\
