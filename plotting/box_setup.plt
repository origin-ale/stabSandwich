# Shared setup for the per-cycle box-plot scripts
# (boxplot_resources_prob_tm_campspp.plt, boxplot_resources_prob_tm_pp.plt).
#
# Set `datafile` before loading this. Optionally set first:
#   maxpanels  panel budget (default 16)
#   logx       1 for a log cycle axis, i.e. log-log panels (default 0)
#
# `logx` is one-shot: box_panels.plt resets it to 0 after drawing, so a log-log
# section cannot leak into the linear sections that follow it.
#
# Mirrors hist_setup.plt: defines blocktitle()/nblocks(), the variables `nb`,
# `i0` and `rows`/`cols`, and sizes the terminal to match that grid, so the
# script adapts to however many p values the run produced.
#
# NOTE: `first` is a reserved coordinate keyword in gnuplot, hence the name `i0`.

if (!exists("maxpanels")) { maxpanels = 16 }
if (!exists("logx")) { logx = 0 }

blocktitle(f, i) = system(sprintf("grep -a '^# p =' '%s' | sed -n '%dp' | sed 's/^# //; s/,.*//'", f, i + 1))

nblocks(f) = system(sprintf("grep -ac '^# p =' '%s'", f)) + 0

nb = nblocks(datafile)

# More p values than the panel budget: drop the lowest ones, as hist_setup does.
i0 = (nb > maxpanels) ? nb - maxpanels : 0
npanels = nb - i0
cols = ceil(sqrt(npanels))
rows = ceil(1.0 * npanels / cols)

eval sprintf("set term png size %d, %d font \",14\"", 380 * cols, 300 * rows)

set style fill solid 0.25 border lc rgb "black"
set grid ytics
set tics font ",11"
set xlabel font ",13"
set ylabel font ",13"
set xlabel "Cycle"

# Resources span orders of magnitude across cycles; a log axis is what makes the
# quartile spread legible at both ends.
set logscale y
set format y "%1g"

if (logx) {
    set logscale x
    set format x "%1g"
    # Auto tics on a log axis crowd into an unreadable smear at the top of the
    # cycle range; gnuplot draws only those of these that fall in range. The
    # list stays dense above 10 because a CAMPS-PP resource that only starts at
    # its switch layer spans well under a decade, and a sparser list would leave
    # such a panel with a single labelled tic.
    set xtics (1, 2, 3, 5, 7, 10, 15, 20, 30, 50)
} else {
    unset logscale x
    set format x "% h"
    set xtics autofreq
}
