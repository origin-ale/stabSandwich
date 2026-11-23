# Box plots of the per-sample CAMPS-PP computation time against the switch point c
# (gates left for PP), one box series per depth M.
#
# Input: M<M>_switch_optimization_full_<N>.txt as written by switch_optimization.jl:
#   c  time_1 … time_Nsamples
# one row per switch point. The box spans the quartiles (Julia's default quantile,
# type 7), the whiskers the extrema over samples, the bar inside is the median and
# the line the mean; the individual samples are drawn behind the boxes.
#
# This is the distribution behind plot_switches.plt, which shows only the mean and
# its standard error (from the _avgs_ file).
#
# dir  input directory (default: the newest fetched run's copy of output/switch_optimization;
#      run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# N    system size (default "20")
# Ms   the M values to plot, space-separated (default: every M with a _full_<N> file in dir)
# e.g. gnuplot -e "run='20261123-091438-switch_optimization'" plotting/plot_switches_box.plt
datadir = "output/switch_optimization"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)
if (!exists("N")) { N = "20" }
if (!exists("Ms")) {
  Ms = system(sprintf("ls '%s' | sed -n 's/^M\\([0-9]*\\)_switch_optimization_full_%s\\.txt$/\\1/p' | sort -n | paste -sd ' ' -", dir, N))
}
nM = words(Ms)
if (nM == 0) {
  print "plot_switches_box.plt: no M*_switch_optimization_full_".N.".txt in ".dir
  exit status 1
}

fullfile(M) = dir.'/M'.M.'_switch_optimization_full_'.N.'.txt'

# Per-row statistics of the samples (columns 2..NF):
#   1 c  2 mean  3 q1  4 median  5 q3  6 min  7 max
# macOS awk has no asort, hence the insertion sort (50 samples per row). The awk
# program is a single-quoted gnuplot string, so no escapes are processed in it.
statawk = '!/^#/ && NF > 1 { n = NF - 1; s = 0; for (i = 1; i <= n; i++) { v = $(i+1); s += v; j = i - 1; while (j > 0 && x[j] > v) { x[j+1] = x[j]; j-- } x[j+1] = v } print $1, s/n, q(x,n,0.25), q(x,n,0.5), q(x,n,0.75), x[1], x[n] } function q(x, n, p,   h, lo) { h = (n-1)*p; lo = int(h); return (lo + 1 < n) ? x[lo+1] + (h-lo)*(x[lo+2]-x[lo+1]) : x[n] }'
statcmd(M) = "< awk -v OFMT=%.9g '".statawk."' '".fullfile(M)."'"

# Number of samples per row (columns 2..nsamp+1), from the first data row.
nsamp(M) = system(sprintf("awk '!/^#/ && NF > 1 { print NF - 1; exit }' '%s'", fullfile(M))) + 0

set term png size 850, 600 font ",16"
if (nM == 1) {
  set output plotdir."/M".Ms."_switch_optimization_box_".N.".png"
  set title "CAMPS-PP computation times, ".N."×(".N."+".Ms.") circuit" font ",24"
} else {
  set output plotdir."/varM_switch_optimization_box_".N.".png"
  set title "CAMPS-PP computation times, ".N."×(".N."+M) circuit" font ",24"
}

set logscale y
set format y "%1g"
set xlabel "c"
set ylabel "Time (s)"
# y range from the sample extrema over every plotted M: log autoscaling rounds out to
# whole decades, which here leaves an empty decade below the fastest samples.
yext = system(sprintf("cat %s | awk '!/^#/ && NF > 1 { for (i = 2; i <= NF; i++) { if (lo == \"\" || $i < lo) lo = $i; if ($i > hi) hi = $i } } END { print lo, hi }'", \
  system("for M in ".Ms."; do printf \"'%s' \" \"".dir."/M${M}_switch_optimization_full_".N.".txt\"; done")))
set yrange [word(yext, 1) / 1.6 : word(yext, 2) * 1.6]
set grid ytics
set key top right
set xtics 1
set style fill solid 0.25 border lc rgb "black"

# Several M: offset the series sideways so the boxes sit side by side.
bw = 0.6 / nM
set boxwidth bw * 0.85 absolute
off(k) = (k - (nM + 1) / 2.0) * bw
set autoscale xfix
set offsets 0.6, 0.6, 0, 0

boxc(k) = word("#4269d0 #efb118 #ff725c #6cc5b0 #3ca951 #a463f2", (k - 1) % 6 + 1)
darkc(k) = word("#1b3a8c #9c7016 #a8402f #3a7d6d #1f6b2e #6a2fb0", (k - 1) % 6 + 1)
# Samples: spread across the box width by sample index, so they don't stack in a column.
jit(k, j, n) = off(k) + (n > 1 ? ((j - 2.0) / (n - 1) - 0.5) * bw * 0.7 : 0)

plot for [k = 1:nM] for [j = 2:nsamp(word(Ms, k)) + 1] fullfile(word(Ms, k)) \
       using ($1 + jit(k, j, nsamp(word(Ms, k)))):j with points pt 7 ps 0.35 \
       lc rgb "#b0b0b0" notitle, \
     for [k = 1:nM] statcmd(word(Ms, k)) using ($1 + off(k)):3:6:7:5 \
       with candlesticks whiskerbars 0.5 lc rgb boxc(k) \
       title (nM == 1 ? "quartiles / extrema" : "M=".word(Ms, k)), \
     for [k = 1:nM] statcmd(word(Ms, k)) using ($1 + off(k)):4:4:4:4 \
       with candlesticks lw 3 lc rgb darkc(k) notitle, \
     for [k = 1:nM] statcmd(word(Ms, k)) using ($1 + off(k)):2 \
       with linespoints lw 2 pt 7 ps 0.8 lc rgb darkc(k) \
       title (nM == 1 ? "mean" : "M=".word(Ms, k)." mean")
