# Box plots of the per-cycle CAMPS-PP resource distributions over samples, one
# panel per magic probability.
#
# Input: output/carlos_varseed64/resources_prob_tm_*_campspp.txt as written by
# save_stats_maxcol. This is the distribution behind
# plot_resources_prob_tm_campspp.plt, which shows only the sample mean and its
# standard error.

# == Number of Pauli strings ============================================================
datafile = "output/carlos_varseed64/resources_prob_tm_NP_campspp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_NP_campspp_box.png"
ylab = "N. of Pauli strings"
plottitle = "46-qubit XXZ circuit — CAMPS-PP Pauli string number"
load "plotting/box_panels.plt"

# == Mean Pauli weight ==================================================================
datafile = "output/carlos_varseed64/resources_prob_tm_PW_campspp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_PW_campspp_box.png"
ylab = "Avg. Pauli weight"
plottitle = "46-qubit XXZ circuit — CAMPS-PP mean Pauli weight"
load "plotting/box_panels.plt"

# == Mean Pauli weight, log-log =========================================================
# The weight grows as a power of the cycle rather than exponentially, so a log
# cycle axis is what straightens it out; cf. plot_resources_prob_tm_campspp.plt,
# which uses `set logscale xy` for this quantity alone.
datafile = "output/carlos_varseed64/resources_prob_tm_PW_campspp.txt"
logx = 1
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_PW_campspp_box_loglog.png"
ylab = "Avg. Pauli weight"
plottitle = "46-qubit XXZ circuit — CAMPS-PP mean Pauli weight (log-log)"
load "plotting/box_panels.plt"

# == Mean |Pauli coefficient| ===========================================================
datafile = "output/carlos_varseed64/resources_prob_tm_PC_campspp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_PC_campspp_box.png"
ylab = "Avg. |Pauli coeff.|"
plottitle = "46-qubit XXZ circuit — CAMPS-PP mean Pauli coeff."
load "plotting/box_panels.plt"

# == Bond dimension =====================================================================
datafile = "output/carlos_varseed64/resources_prob_tm_bd_campspp.txt"
load "plotting/box_setup.plt"
set output "output/Resources_prob_tm_bd_campspp_box.png"
ylab = "Bond dimension"
plottitle = "46-qubit XXZ circuit — CAMPS-PP bond dimension"
load "plotting/box_panels.plt"
