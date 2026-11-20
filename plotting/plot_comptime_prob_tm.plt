# Mean computation time (with its standard error) against magic probability,
# CAMPS-PP against PP.
#
# Input: comptime_prob_tm.txt as written by Comptime_prob_tm.jl, with seven
# columns per method (1 p, then per method mean, std.err, q1, median, q3, min,
# max), so method j (0-based, in the `methods` order of the run) starts at
# column 2 + 7*j.
#
# dir_campspp, dir_pp  directories holding each method's comptime_prob_tm.txt
#                      (default: output/Comptime_prob_tm for both, i.e. one run
#                      with methods=campspp,pp)
# j_campspp, j_pp      each method's 0-based position in its file's `methods`
#                      (default 0 and 1)
# Old data (one method per run):
#   gnuplot -e "dir_campspp='output/old/carlos_varseed64'; dir_pp='output/old/carlos_varseed_pp'; j_pp=0" plotting/plot_comptime_prob_tm.plt
if (!exists("dir_campspp")) { dir_campspp = "output/Comptime_prob_tm" }
if (!exists("dir_pp")) { dir_pp = "output/Comptime_prob_tm" }
if (!exists("j_campspp")) { j_campspp = 0 }
if (!exists("j_pp")) { j_pp = 1 }
col(j, q) = 2 + 7*j + q

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

plot dir_campspp.'/comptime_prob_tm.txt' using 1:col(j_campspp,0):col(j_campspp,1) with yerrorlines title "CAMPS-PP, χ=64, thl=10^{-10}",\
dir_pp.'/comptime_prob_tm.txt' using 1:col(j_pp,0):col(j_pp,1) with yerrorlines title "PP, thl 10^{-10}"
# Old χ = 128 and χ = 32 runs: output/old/carlos_varseed and output/old/carlos_varseed32
# input using 1:4:5 with yerrorlines title "CAMPS - search"\,
# input using 1:6:7 with yerrorlines title "MPS (thl 10^{-9})",\
