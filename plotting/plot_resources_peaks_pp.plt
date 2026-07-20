# Histograms of the per-sample peak resource values (one panel per magic probability p).
# Input: output/carlos_varseed_pp/resources_prob_tm_*_peak_pp.txt
#        column 1 = sample index, column 2 = peak value; one data block per p.
#
# The panel grid and the x-ranges follow the data, so runs with more p values plot
# without edits. If there are more p values than `maxpanels`, the lowest ones are
# dropped — raise it here, or override with `gnuplot -e "maxpanels=N"`.
if (!exists("maxpanels")) { maxpanels = 16 }

# == Number of Pauli strings ============================================================
# Peaks span ~47 → ~7e5, so histogram the base-10 logarithm.
datafile = "output/carlos_varseed_pp/resources_prob_tm_NP_peak_pp.txt"
load "plotting/hist_setup.plt"
set output "output/Peaks_prob_tm_NP_pp_hist.png"
set xlabel "log_{10}(peak n. of Pauli strings)"
binexpr = "log10($2)"
w = 0.25
plottitle = sprintf("46-qubit XXZ circuit — PP peak Pauli string number (%d samples per p)", nsamples(datafile))
load "plotting/hist_panels.plt"

# == Mean Pauli weight ==================================================================
datafile = "output/carlos_varseed_pp/resources_prob_tm_PW_peak_pp.txt"
load "plotting/hist_setup.plt"
set output "output/Peaks_prob_tm_PW_pp_hist.png"
set xlabel "Peak avg. Pauli weight"
binexpr = "$2"
w = 0.5
plottitle = sprintf("46-qubit XXZ circuit — PP peak mean Pauli weight (%d samples per p)", nsamples(datafile))
load "plotting/hist_panels.plt"

# == Mean |Pauli coefficient| ===========================================================
# Values are multiples of 1/47; bin on that natural spacing.
datafile = "output/carlos_varseed_pp/resources_prob_tm_PC_peak_pp.txt"
load "plotting/hist_setup.plt"
set output "output/Peaks_prob_tm_PC_pp_hist.png"
set xlabel "Peak avg. |Pauli coeff.|"
binexpr = "$2"
w = 1.0 / 47.0
plottitle = sprintf("46-qubit XXZ circuit — PP peak mean |Pauli coeff.| (%d samples per p)", nsamples(datafile))
load "plotting/hist_panels.plt"
