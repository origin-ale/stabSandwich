# Draws one multiplot page: a per-cycle box plot panel per plotted block.
#
# Expects box_setup.plt to have been loaded (for `nb`, `i0`, `rows`, `cols`),
# plus `plottitle` and `ylab`.
#
# Column layout written by save_stats_maxcol:
#   1 layer  2 mean  3 std.err  4 max-sample  5 n-samples
#   6 q1     7 median  8 q3     9 min        10 max
# The box spans the quartiles, the whiskers the extrema over samples, and the
# mean is overlaid so it can be read against the median — the gap between the
# two is the asymmetry that a mean-and-error-bar plot hides.

set ylabel ylab

# One y-range spanning every plotted block, so the panels stay comparable.
# Range over the whiskers (min and max), not the boxes, or the extrema clip.
set autoscale
eval sprintf("stats \"%s\" index %d:%d using 9 nooutput", datafile, i0, nb - 1)
ymin = STATS_min
eval sprintf("stats \"%s\" index %d:%d using 10 nooutput", datafile, i0, nb - 1)
set yrange [ymin * 0.8 : STATS_max * 1.25]

eval sprintf("set multiplot layout %d,%d title \"%s\" font \",18\"", rows, cols, plottitle)
do for [i = i0:nb-1] {
    set title blocktitle(datafile, i) font ",14"
    # Resources rise steeply with the cycle, so no corner is reliably free
    # across every p; `opaque` keeps the key readable wherever the curve runs.
    set key font ",10" top left opaque box
    plot datafile index i using 1:6:9:10:8 with candlesticks whiskerbars 0.5 \
           lc rgb "#4269d0" title "quartiles / extrema", \
         "" index i using 1:7:7:7:7 with candlesticks lw 3 \
           lc rgb "#1b3a8c" notitle, \
         "" index i using 1:2 with linespoints lw 2 pt 7 ps 0.5 \
           lc rgb "#e15759" title "mean"
}
unset multiplot
