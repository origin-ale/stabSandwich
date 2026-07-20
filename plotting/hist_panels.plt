# Draws one multiplot page: a histogram panel per plotted block.
#
# Expects hist_setup.plt to have been loaded (for `nb`, `i0`, `rows`, `cols`),
# plus `binexpr` — the quantity to histogram, written as a using-expression
# string — `w` (bin width) and `plottitle`.

# One x-range spanning every plotted block, so the panels stay comparable.
# `stats` silently ignores samples outside the *active* ranges, and the count
# axis set above is [0:nsmp] — which would throw away every negative log10 value
# — so autoscale across the measurement, then restore the count axis.
set autoscale
eval sprintf("stats \"%s\" index %d:%d using (%s) nooutput", datafile, i0, nb - 1, binexpr)
set yrange [0:nsmp]
set xrange [bin(STATS_min, w) - w : bin(STATS_max, w) + w]

# Absolute box width: a relative one collapses to nothing when a block's values
# all land in a single bin (as the peak counts do at p = 0).
set boxwidth w

eval sprintf("set multiplot layout %d,%d title \"%s\" font \",18\"", rows, cols, plottitle)
do for [i = i0:nb-1] {
    set title blocktitle(datafile, i) font ",14"
    eval sprintf("plot datafile index %d using (bin(%s, %.10g)):(1) smooth freq with boxes lc %d", i, binexpr, w, i + 1)
}
unset multiplot
