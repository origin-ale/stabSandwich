# Box plots of the per-sample computation time against t/N, exact CAMPS against
# exact PP, one box series per method side by side.
#
# Input: comptimes_ratios_<N>_full.txt as written by Comptime_ratios.jl:
#   t/N  time_1 … time_Nsamples
# one row per ratio, one blank-line-separated block per method in the order of the
# header's methods=, then one block for the gate conversion (a sub-timing of CAMPS).
# The box spans the quartiles (Julia's default quantile, type 7), the whiskers the
# extrema over samples, the bar inside is the median and the line the mean; the
# individual samples are drawn behind the boxes.
#
# Comptime_ratios.jl records no per-sample compile time, so nothing is subtracted:
# a sample where the warm-up missed keeps its compilation (the run's log has a
# WARNING per such sample; their count is printed).
#
# This is the distribution behind plot_ratios_times.plt, which shows only the mean
# and its standard error (from the _avgs_ file).
#
# dir      input directory (default: the newest fetched run's copy of output/Comptime_ratios;
#          run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# N        system size (default: the one in dir's _full_ files; the largest if there are several)
# methods  the blocks to plot, space-separated: any of the run's methods, and convert
#          (default: the run's methods, as in the file header)
# e.g. gnuplot -e "run='20261125-163541-Comptime_ratios'; methods='camps pp convert'" plotting/plot_ratios_box.plt
datadir = "output/Comptime_ratios"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)
if (!exists("N")) {
  Ns = system(sprintf("ls '%s' | sed -n 's/^comptimes_ratios_\\([0-9]*\\)_full\\.txt$/\\1/p' | sort -nu | paste -sd ' ' -", dir))
  if (words(Ns) == 0) {
    print "plot_ratios_box.plt: no comptimes_ratios_*_full.txt in ".dir
    exit status 1
  }
  N = word(Ns, words(Ns))
  if (words(Ns) > 1) { print "plot_ratios_box.plt: ".dir." has N = ".Ns."; plotting N = ".N." (pass N=… for another)" }
}
full = dir.'/comptimes_ratios_'.N.'_full.txt'
if (system(sprintf("[ -f '%s' ] && echo 1 || echo 0", full)) + 0 == 0) {
  print "plot_ratios_box.plt: no ".full
  exit status 1
}

# The run's methods, in block order (methods=camps,pp, or Julia's methods=[:camps, :pp])
runmethods = system(sprintf('head -1 "%s" | grep -oE "methods=(\[[^]]*\]|[^ ]*)" | cut -d= -f2 | tr -d "[]: " | tr , " "', full))
if (words(runmethods) == 0) {
  print "plot_ratios_box.plt: no methods= in the header of ".full
  exit status 1
}
if (!exists("methods")) { methods = runmethods }
nM = words(methods)

# Data rows of block b (0-based; blocks are separated by single blank lines)
blockawk = 'NF == 0 { if (seen) { blk++; seen = 0 } next } /^#/ { next } { seen = 1; if (blk == b) print }'
nblocks = system("awk 'NF == 0 { if (seen) { n++; seen = 0 } next } /^#/ { next } { seen = 1 } END { print n + seen }' '".full."'") + 0
# convert is the block after the methods'
block(m) = (m eq "convert") ? words(runmethods) : method_pos(full, m)
do for [m in methods] {
  if (block(m) < 0 || block(m) >= nblocks) {
    print "plot_ratios_box.plt: no ".m." block in ".full." (methods: ".runmethods.", convert)"
    exit status 1
  }
}
samplescmd(m) = sprintf("awk -v b=%d '%s' '%s'", block(m), blockawk, full)

# Per-row statistics of the samples (columns 2..NF):
#   1 t/N  2 mean  3 q1  4 median  5 q3  6 min  7 max
# macOS awk has no asort, hence the insertion sort (100 samples per row). The awk
# program is a single-quoted gnuplot string, so no escapes are processed in it.
statawk = '{ n = NF - 1; s = 0; for (i = 1; i <= n; i++) { v = $(i+1); s += v; j = i - 1; while (j > 0 && x[j] > v) { x[j+1] = x[j]; j-- } x[j+1] = v } print $1, s/n, q(x,n,0.25), q(x,n,0.5), q(x,n,0.75), x[1], x[n] } function q(x, n, p,   h, lo) { h = (n-1)*p; lo = int(h); return (lo + 1 < n) ? x[lo+1] + (h-lo)*(x[lo+2]-x[lo+1]) : x[n] }'
statcmd(m) = "< ".samplescmd(m)." | awk -v OFMT=%.9g '".statawk."'"

# Ratio spacing, from the first two rows: sets the box widths
dx = system(samplescmd(word(methods, 1))." | awk 'NR == 1 { a = $1 } NR == 2 { print $1 - a; exit }'") + 0
if (dx <= 0) { dx = 0.05 }

# Warm-up misses, from the run's log (output/runs/<id>/log.txt) when there is one
logfile = system(sprintf("echo '%s' | sed -n 's|^\\(output/runs/[^/]*\\)/.*|\\1/log.txt|p'", dir))
if (strlen(logfile) > 0) {
  nmiss = system(sprintf("grep -c 'warm-up missed' '%s' 2>/dev/null || true", logfile)) + 0
  if (nmiss > 0) {
    print sprintf("plot_ratios_box.plt: %d timed samples include compilation (warm-up missed, see %s)", nmiss, logfile)
  }
}

label(m) = (m eq "camps") ? "Exact CAMPS" : (m eq "pp") ? "Exact PP" : (m eq "convert") ? "Gate conversion" : m

set term png size 850, 600 font ",16"
set output plotdir."/CompTimes_ratios_".N."_box.png"
set title "Computation times, N=".N." t-rotation circuit" font ",24"

set logscale y
set format y "%1g"
set xlabel "t/N"
set ylabel "Time (s)"
# y range from the sample extrema over every plotted method: log autoscaling rounds
# out to whole decades, which can leave an empty decade below the fastest samples
yext = ""
do for [m in methods] { yext = yext." ".system(samplescmd(m)." | awk '{ for (i = 2; i <= NF; i++) if ($i > 0) { if (lo == \"\" || $i < lo) lo = $i; if ($i > hi) hi = $i } } END { print lo, hi }'") }
ylo = 0; yhi = 0
do for [i = 1:words(yext):2] {
  if (ylo == 0 || word(yext, i) + 0 < ylo) { ylo = word(yext, i) + 0 }
  if (word(yext, i + 1) + 0 > yhi) { yhi = word(yext, i + 1) + 0 }
}
set yrange [ylo / 1.6 : yhi * 1.6]
set grid ytics
set key top left
set xtics dx format "%.2f"
set style fill solid 0.25 border lc rgb "black"

# The methods' boxes sit side by side around each ratio
bw = 0.75 * dx / nM
set boxwidth bw * 0.85 absolute
off(k) = (k - (nM + 1) / 2.0) * bw
set autoscale xfix
set offsets 0.6 * dx, 0.6 * dx, 0, 0

boxc(k) = word("#4269d0 #ff725c #6cc5b0 #efb118", (k - 1) % 4 + 1)
darkc(k) = word("#1b3a8c #a8402f #3a7d6d #9c7016", (k - 1) % 4 + 1)
# Samples, spread across the box width by sample index so they don't stack in a column
jitawk = '{ n = NF - 1; for (i = 2; i <= NF; i++) print $1 + o + (n > 1 ? ((i - 2) / (n - 1) - 0.5) * w : 0), $i }'
jitcmd(k) = "< ".samplescmd(word(methods, k))." | awk -v o=".sprintf("%.9g", off(k))." -v w=".sprintf("%.9g", bw * 0.7)." '".jitawk."'"

plot for [k = 1:nM] jitcmd(k) using 1:2 with points pt 7 ps 0.3 lc rgb "#b0b0b0" notitle, \
     for [k = 1:nM] statcmd(word(methods, k)) using ($1 + off(k)):3:6:7:5 \
       with candlesticks whiskerbars 0.5 lc rgb boxc(k) title label(word(methods, k)), \
     for [k = 1:nM] statcmd(word(methods, k)) using ($1 + off(k)):4:4:4:4 \
       with candlesticks lw 3 lc rgb darkc(k) notitle, \
     for [k = 1:nM] statcmd(word(methods, k)) using ($1 + off(k)):2 \
       with linespoints lw 2 pt 7 ps 0.8 lc rgb darkc(k) title label(word(methods, k))." mean"
