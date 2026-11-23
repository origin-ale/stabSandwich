# Box plots of the per-cycle CAMPS bond-dimension distribution over samples, one
# panel per doping count n (or magic probability p, in p-based data). The bond
# dimension is the only resource CAMPS records.
#
# Input: resources_prob_tm_bd_camps.txt as written by save_stats_maxcol
# (the CAMPS-PP counterpart is the last section of
# boxplot_resources_prob_tm_campspp.plt).
#
# Input directory `dir` (default: the newest fetched run's copy of output/Comptime_prob_tm;
# run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# Old data: gnuplot -e "dir='output/old/carlos_indep_camps'" plotting/boxplot_resources_prob_tm_camps.plt
datadir = "output/Comptime_prob_tm"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)

# == Bond dimension =====================================================================
datafile = dir."/resources_prob_tm_bd_camps.txt"
load "plotting/box_setup.plt"
set output plotdir."/Resources_prob_tm_bd_camps_box.png"
ylab = "Bond dimension"
plottitle = "46-qubit XXZ circuit — CAMPS bond dimension"
load "plotting/box_panels.plt"
