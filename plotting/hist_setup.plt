# Shared setup for the per-sample histogram scripts (plot_resources_peaks_pp.plt,
# plot_comptime_peaks_pp.plt).
#
# Set `datafile` before loading this. Optionally set `maxpanels` first to change
# the panel budget (default 16).
#
# Defines the helpers bin()/blocktitle(), and the variables `nb` (number of data
# blocks in the file), `i0` (index of the first block to plot) and `rows`/`cols`
# (the multiplot grid). Sizes the terminal to match that grid, so the script
# adapts to however many p values the run produced.
#
# NOTE: `first` is a reserved coordinate keyword in gnuplot, hence the name `i0`.

if (!exists("maxpanels")) { maxpanels = 16 }

# Title for data index i: pulls the matching "# p = …, μ = …" header line
# from the file (data index i corresponds to the (i+1)-th such header line),
# keeping only the "p = …" field.
blocktitle(f, i) = system(sprintf("grep -a '^# p =' '%s' | sed -n '%dp' | sed 's/^# //; s/,.*//'", f, i + 1))

nblocks(f) = system(sprintf("grep -ac '^# p =' '%s'", f)) + 0

# Samples per block, read off the "samples=…" field of the run header.
nsamples(f) = system(sprintf("grep -a -m1 -o 'samples=[0-9]*' '%s' | cut -d= -f2", f)) + 0

# Bin centre for value x with bin width w.
bin(x, w) = w * (floor(x / w) + 0.5)

nb = nblocks(datafile)

# More p values than the panel budget: drop the lowest ones. Those are the least
# informative — at small p every sample piles into a single bin — and the blocks
# are ordered by increasing p, so dropping a prefix keeps the interesting tail.
i0 = (nb > maxpanels) ? nb - maxpanels : 0
npanels = nb - i0
cols = ceil(sqrt(npanels))
rows = ceil(1.0 * npanels / cols)

eval sprintf("set term png size %d, %d font \",14\"", 380 * cols, 300 * rows)

set style fill solid 0.55 border lc rgb "black"
set grid ytics
unset key
set tics font ",11"
set xlabel font ",13"
set ylabel font ",13"
set ylabel "N. of samples"

# Common count axis: every block holds the same number of samples, so a shared
# y-range is what makes the panels comparable across p.
nsmp = nsamples(datafile)
set yrange [0:nsmp]
