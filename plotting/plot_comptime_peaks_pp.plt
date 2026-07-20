# Histograms of the per-sample computation times (one panel per magic probability p).
# Input: output/carlos_varseed_pp/comptime_full_prob_tm_pp.txt
#        column 1 = sample index, column 2 = wall time (s), column 3 = GC time (s);
#        one data block per p.
#
# The panel grid and the x-ranges follow the data, so runs with more p values plot
# without edits. If there are more p values than `maxpanels`, the lowest ones are
# dropped — raise it here, or override with `gnuplot -e "maxpanels=N"`.
if (!exists("maxpanels")) { maxpanels = 16 }

datafile = "output/carlos_varseed_pp/comptime_full_prob_tm_pp.txt"
load "plotting/hist_setup.plt"

# == Total wall time ====================================================================
# Times span ~0.07 s → ~17 s, so histogram the base-10 logarithm.
set output "output/Peaks_prob_tm_time_pp_hist.png"
set xlabel "log_{10}(computation time / s)"
binexpr = "log10($2)"
w = 0.2
plottitle = sprintf("46-qubit XXZ circuit — PP computation time (%d samples per p)", nsamples(datafile))
load "plotting/hist_panels.plt"

# == GC time as a fraction of the total =================================================
# Most samples do no GC at all, so this is plotted linearly: the bin at 0 counts
# those runs, and a log axis could not represent them.
set output "output/Peaks_prob_tm_gcfrac_pp_hist.png"
set xlabel "GC time (% of total)"
binexpr = "100.0 * $3 / $2"
w = 2.5
plottitle = sprintf("46-qubit XXZ circuit — PP garbage-collection share (%d samples per p)", nsamples(datafile))
load "plotting/hist_panels.plt"
