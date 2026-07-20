# Draws one multiplot page: a per-cycle box plot panel per plotted block.
#
# Expects box_setup.plt to have been loaded (for `nb`, `i0`, `rows`, `cols`,
# `logx`), plus `plottitle` and `ylab`.
#
# Column layout written by save_stats_maxcol:
#   1 layer  2 mean  3 std.err  4 max-sample  5 n-samples
#   6 q1     7 median  8 q3     9 min        10 max
# The box spans the quartiles, the whiskers the extrema over samples, and the
# mean is overlaid so it can be read against the median — the gap between the
# two is the asymmetry that a mean-and-error-bar plot hides.

set ylabel ylab

# Ranges spanning every plotted block, so the panels stay comparable.
#
# All three measurements are taken while fully autoscaled and only then applied:
# `stats` silently drops records lying outside the *active* ranges, so setting
# the y-range before measuring the cycles truncates the measured cycle span.
set autoscale
eval sprintf("stats \"%s\" index %d:%d using 9 nooutput", datafile, i0, nb - 1)
ymin = STATS_min
eval sprintf("stats \"%s\" index %d:%d using 10 nooutput", datafile, i0, nb - 1)
ymax = STATS_max
eval sprintf("stats \"%s\" index %d:%d using 1 nooutput", datafile, i0, nb - 1)
cyclemin = STATS_min
cyclemax = STATS_max

# Range over the whiskers (min and max), not the boxes, or the extrema clip.
set yrange [ymin * 0.8 : ymax * 1.25]

# Cycle 0 is a valid record (PP tracks from the first layer) but has no place on
# a log axis: gnuplot drops it, and clamping the lower limit to 1 only keeps
# that from dragging the axis to zero. Log-log panels therefore start at cycle 1.
# The limits go through named variables because `:` inside `set xrange [a:b]`
# is the range separator, so a ternary written in there is misparsed.
xlo = (logx && cyclemin < 1) ? 1 : cyclemin
xmin = logx ? xlo * 0.9 : xlo - 1
xmax = logx ? cyclemax * 1.15 : cyclemax + 1
set xrange [xmin:xmax]

# Box width, supplied per record as a sixth column. A constant width in x would
# look increasingly pinched towards the right of a log axis, so on log x the
# width is made proportional to x, which is constant in log space.
wexpr = logx ? '($1*0.055)' : '(0.6)'

eval sprintf("set multiplot layout %d,%d title \"%s\" font \",18\"", rows, cols, plottitle)
do for [i = i0:nb-1] {
    set title blocktitle(datafile, i) font ",14"
    # Resources rise steeply with the cycle, so no corner is reliably free
    # across every p; `opaque` keeps the key readable wherever the curve runs.
    set key font ",10" top left opaque box
    eval sprintf('plot datafile index %d using 1:6:9:10:8:%s \
                    with candlesticks whiskerbars 0.5 \
                    lc rgb "#4269d0" title "quartiles / extrema", \
                  "" index %d using 1:7:7:7:7:%s with candlesticks lw 3 \
                    lc rgb "#1b3a8c" notitle, \
                  "" index %d using 1:2 with linespoints lw 2 pt 7 ps 0.5 \
                    lc rgb "#e15759" title "mean"', i, wexpr, i, wexpr, i)
}
unset multiplot

# `logx` is one-shot (see box_setup.plt): reset so the next section starts linear.
logx = 0
