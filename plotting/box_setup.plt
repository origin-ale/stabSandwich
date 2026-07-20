# Shared setup for the per-cycle box-plot scripts (plot_resources_prob_tm_box.plt).
#
# Set `datafile` before loading this. Optionally set `maxpanels` first to change
# the panel budget (default 16).
#
# Mirrors hist_setup.plt: defines blocktitle()/nblocks(), the variables `nb`,
# `i0` and `rows`/`cols`, and sizes the terminal to match that grid, so the
# script adapts to however many p values the run produced.
#
# NOTE: `first` is a reserved coordinate keyword in gnuplot, hence the name `i0`.

if (!exists("maxpanels")) { maxpanels = 16 }

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

# Boxes sit on integer cycles, so an absolute width just under 1 keeps
# neighbouring boxes from touching at any panel size.
set boxwidth 0.6 absolute

# Resources span orders of magnitude across cycles; a log axis is what makes the
# quartile spread legible at both ends.
set logscale y
set format y "%1g"
