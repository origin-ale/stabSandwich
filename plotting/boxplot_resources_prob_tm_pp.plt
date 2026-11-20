# Box plots of the per-cycle PP resource distributions over samples, one panel
# per magic probability.
#
# Input: resources_prob_tm_*_pp.txt as written by
# save_stats_maxcol. This is the distribution behind
# plot_resources_prob_tm_pp.plt, which shows only the sample mean and its
# standard error.
#
# Input directory `dir` (default: output/Comptime_prob_tm, where Comptime_prob_tm.jl
# writes). Old data: gnuplot -e "dir='output/old/carlos_varseed_pp'" plotting/boxplot_resources_prob_tm_pp.plt
if (!exists("dir")) { dir = "output/Comptime_prob_tm" }

# == Number of Pauli strings ============================================================
datafile = dir."/resources_prob_tm_NP_pp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_NP_pp_box.png"
ylab = "N. of Pauli strings"
plottitle = "46-qubit XXZ circuit — PP Pauli string number"
load "plotting/box_panels.plt"

# == Mean Pauli weight ==================================================================
datafile = dir."/resources_prob_tm_PW_pp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_PW_pp_box.png"
ylab = "Avg. Pauli weight"
plottitle = "46-qubit XXZ circuit — PP mean Pauli weight"
load "plotting/box_panels.plt"

# == Mean Pauli weight, log-log =========================================================
# The weight grows as a power of the cycle rather than exponentially, so a log
# cycle axis is what straightens it out; cf. plot_resources_prob_tm_pp.plt,
# which uses `set logscale xy` for this quantity alone.
datafile = dir."/resources_prob_tm_PW_pp.txt"
logx = 1
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_PW_pp_box_loglog.png"
ylab = "Avg. Pauli weight"
plottitle = "46-qubit XXZ circuit — PP mean Pauli weight (log-log)"
load "plotting/box_panels.plt"

# == Mean |Pauli coefficient| ===========================================================
datafile = dir."/resources_prob_tm_PC_pp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_PC_pp_box.png"
ylab = "Avg. |Pauli coeff.|"
plottitle = "46-qubit XXZ circuit — PP mean Pauli coeff."
load "plotting/box_panels.plt"
