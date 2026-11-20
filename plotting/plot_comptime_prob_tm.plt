# Mean computation time (with its standard error) against the number of doped
# triples per sample (or the magic probability, in p-based data),
# CAMPS-PP against PP.
#
# Input: comptime_prob_tm.txt as written by Comptime_prob_tm.jl, with seven
# columns per method (1 n or p, then per method mean, std.err, q1, median, q3, min,
# max), so method j (0-based, in the `methods` order of the run) starts at
# column 2 + 7*j.
#
# dir_campspp, dir_pp  directories holding each method's comptime_prob_tm.txt
#                      (default for both: the newest fetched run's copy of
#                      output/Comptime_prob_tm, i.e. one run with methods=campspp,pp;
#                      run=<id> picks another, see run.plt)
# Plots go to output/plots/<run id> (<id1>+<id2> for dirs from two runs).
# j_campspp, j_pp      each method's 0-based position in its file's `methods`
#                      (default 0 and 1)
# Old data (one method per run):
#   gnuplot -e "dir_campspp='output/old/carlos_varseed64'; dir_pp='output/old/carlos_varseed_pp'; j_pp=0" plotting/plot_comptime_prob_tm.plt
datadir = "output/Comptime_prob_tm"
load "plotting/run.plt"
if (!exists("dir_campspp")) { dir_campspp = rundir }
if (!exists("dir_pp")) { dir_pp = rundir }
plotdir = plot_dir(dir_campspp." ".dir_pp)
if (!exists("j_campspp")) { j_campspp = 0 }
if (!exists("j_pp")) { j_pp = 1 }
col(j, q) = 2 + 7*j + q
# x axis: n (doped triples per sample) for count-mode runs, p for p-based data
xlab(f) = (system(sprintf("grep -ac 'magic_counts=' '%s'", f)) + 0 > 0) ? "Doped triples per sample" : "Magic probability"

set term png size 850, 600 font ",16"
set output plotdir."/CompTimes_prob_tm_graph.png"

set linetype 3 lc rgb "dark-red"

set title "46-qubit 23-layer XXZ circuit" font ",24"
set logscale y
set format y "%1g"
set xlabel xlab(dir_campspp.'/comptime_prob_tm.txt')
set ylabel "Time (s)"

set xrange [0:*]
# Lower end autoscaled: a 0 on a log axis is silently replaced by 0.1, which
# hid PP's sub-0.1 s times at low n.
set yrange [*:1000]

set key top left

plot dir_campspp.'/comptime_prob_tm.txt' using 1:col(j_campspp,0):col(j_campspp,1) with yerrorlines title "CAMPS-PP, χ=64, thl=10^{-10}",\
dir_pp.'/comptime_prob_tm.txt' using 1:col(j_pp,0):col(j_pp,1) with yerrorlines title "PP, thl 10^{-10}"
# Old χ = 128 and χ = 32 runs: output/old/carlos_varseed and output/old/carlos_varseed32
# input using 1:4:5 with yerrorlines title "CAMPS - search"\,
# input using 1:6:7 with yerrorlines title "MPS (thl 10^{-9})",\
