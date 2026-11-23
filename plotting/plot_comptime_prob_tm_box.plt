# Box plots of the per-sample computation time against the number of doped
# triples per sample (or the magic probability, in p-based data).
#
# Input: comptime_prob_tm.txt (and gctime_prob_tm.txt) as written by
# Comptime_prob_tm.jl, with seven columns per method:
#   1 n or p, then per method  mean  std.err  q1  median  q3  min  max
# so method j (0-based, in the `methods` order of the run) starts at column
# 2 + 7*j. The box spans the quartiles, the whiskers the extrema over samples.
#
# This is the distribution behind plot_comptime_prob_tm.plt, which shows only
# the mean and its standard error.
#
# Methods: CAMPS-PP, CAMPS and PP, side by side at each n; a method the run
# didn't include is left out.
# dir_campspp, dir_camps, dir_pp  directories holding each method's files
#     (default for all: the newest fetched run's copy of output/Comptime_prob_tm;
#     run=<id> picks another, see run.plt). Once any is given, only the methods
#     given one are plotted.
# j_campspp, j_camps, j_pp  each method's 0-based position in its file's
#     `methods`; -1 leaves the method out (default: read from the file's
#     methods= header). See comptime_methods.plt.
# Plots go to output/plots/<run id> (<id1>+<id2> for dirs from two runs).
# Old data (one method per run):
#   gnuplot -e "dir_campspp='output/old/carlos_indep'; dir_camps='output/old/carlos_indep_camps'; dir_pp='output/old/carlos_indep_pp'" plotting/plot_comptime_prob_tm_box.plt
datadir = "output/Comptime_prob_tm"
load "plotting/run.plt"
load "plotting/comptime_methods.plt"
# x axis: n (doped triples per sample) for count-mode runs, p for p-based data
xlab(f) = (system(sprintf("grep -ac 'magic_counts=' '%s'", f)) + 0 > 0) ? "Doped triples per sample" : "Magic probability"

set term png size 850, 600 font ",16"

set title "46-qubit 23-layer XXZ circuit" font ",24"
set logscale y
set format y "%1g"
set ylabel "Time (s)"

set xrange [*:*]
set yrange [*:*]
set key top left

set style fill solid 0.25 border lc rgb "black"

# Column of quantity `q` for method j (0-based) in a file with 7 columns/method.
col(j, q) = 2 + 7*j + q
Q1 = 2; MED = 3; Q3 = 4; MIN = 5; MAX = 6; MEAN = 0

# Per-method helpers: file `f` ("comptime" or "gctime") of method m, its
# colours, label and truncation threshold (jof(m), its column offset, comes
# from comptime_methods.plt).
fileof(m, f) = dirof(m).'/'.f.'_prob_tm.txt'
boxc(m) = pick(m, "#4269d0", "#3ca951", "#efb118")
darkc(m) = pick(m, "#1b3a8c", "#1d6b30", "#9c7016")
label(m) = pick(m, "CAMPS-PP", "CAMPS", "PP")
thlof(m) = pick(m, "10^{-10}", "10^{-12}", "10^{-10}")
first_file = fileof(word(mlist, 1), "comptime")
set xlabel xlab(first_file)

# The p spacing is 0.005, so boxes must be far narrower than that; `dx` also
# offsets the two methods sideways so their boxes sit side by side instead of
# on top of each other. For count-mode runs both scale with the smallest gap
# between n values (the same fractions: 0.24 and 0.18 of it).
bw = 0.0012
dx = 0.0009
if (system(sprintf("grep -ac 'magic_counts=' '%s'", first_file)) + 0 > 0) {
  xstep = system(sprintf("awk '!/^#/ && NF {x[k++]=$1} END {m=0; for (i=1;i<k;i++) {d=x[i]-x[i-1]; if (d>0 && (m==0 || d<m)) m=d}; print (m>0 ? m : 1)}' '%s'", first_file)) + 0
  bw = 0.24*xstep
  dx = 0.18*xstep
}
set boxwidth bw absolute
# Sideways offset of each method, centred on the tick: none for a lone method,
# ±dx for two, and for three a spacing of 1.25 box widths (two methods'
# 2*dx = 1.5 box widths would push the outer boxes into the next n's). A lone
# method's end boxes sit on the autoscaled x edges, cut in half, so pad the
# range by a box width (offset boxes already clear them).
nm = words(mlist)
sp = (nm == 3) ? 1.25*bw : 2*dx
do for [i = 1:nm] {
  eval sprintf("off_%s = %.17g", word(mlist, i), (nm < 2) ? 0 : (i - (nm + 1)/2.0)*sp)
}
off(m) = value("off_".m)
if (nm < 2) { set autoscale xfix; set offsets bw, bw, 0, 0 }

# == Wall-clock time ====================================================================
set output plotdir."/CompTimes_prob_tm_box.png"

plot for [m in mlist] fileof(m, "comptime") \
       using ($1+off(m)):col(jof(m),Q1):col(jof(m),MIN):col(jof(m),MAX):col(jof(m),Q3) \
       with candlesticks whiskerbars 0.5 lc rgb boxc(m) title label(m).", thl ".thlof(m), \
     for [m in mlist] fileof(m, "comptime") \
       using ($1+off(m)):col(jof(m),MED):col(jof(m),MED):col(jof(m),MED):col(jof(m),MED) \
       with candlesticks lw 3 lc rgb darkc(m) notitle, \
     for [m in mlist] fileof(m, "comptime") using ($1+off(m)):col(jof(m),MEAN) \
       with linespoints pt 7 ps 0.8 lc rgb darkc(m) title label(m)." mean"

# == GC time ============================================================================
# PP's GC time is 0 in most samples, which a log axis can't show: such boxes
# are drawn clipped at the bottom, and autoscaling ignores them. If no plotted
# box has a positive minimum (e.g. a PP-only run), autoscaling has nothing to
# go on, so take the range from the positive values; if there are none, skip.
# posrange(m) = "<any box with min > 0> <smallest positive> <largest>" over
# method m's q1 … max columns.
posrange(m) = system(sprintf("awk -v a=%d -v b=%d -v mn=%d '!/^#/ && NF { if ($mn > 0) ok = 1; for (i = a; i <= b; i++) if ($i > 0) { if (lo == \"\" || $i < lo) lo = $i; if ($i > hi) hi = $i } } END { print ok + 0, lo + 0, hi + 0 }' '%s'", col(jof(m),Q1), col(jof(m),MAX), col(jof(m),MIN), fileof(m, "gctime")))
gc_ok = 0; gc_lo = 0; gc_hi = 0
do for [m in mlist] {
  r = posrange(m)
  if (word(r, 1) + 0) { gc_ok = 1 }
  if (word(r, 2) + 0 > 0 && (gc_lo == 0 || word(r, 2) + 0 < gc_lo)) { gc_lo = word(r, 2) + 0 }
  if (word(r, 3) + 0 > gc_hi) { gc_hi = word(r, 3) + 0 }
}

if (gc_hi == 0) {
  print "plot_comptime_prob_tm_box.plt: every GC time is 0, nothing to show on a log axis; GC plot skipped"
} else {
  if (!gc_ok) { set yrange [gc_lo/2 : gc_hi*2] }
  set output plotdir."/GCTimes_prob_tm_box.png"
  set ylabel "GC time (s)"

  plot for [m in mlist] fileof(m, "gctime") \
         using ($1+off(m)):col(jof(m),Q1):col(jof(m),MIN):col(jof(m),MAX):col(jof(m),Q3) \
         with candlesticks whiskerbars 0.5 lc rgb boxc(m) title label(m).", thl ".thlof(m), \
       for [m in mlist] fileof(m, "gctime") \
         using ($1+off(m)):col(jof(m),MED):col(jof(m),MED):col(jof(m),MED):col(jof(m),MED) \
         with candlesticks lw 3 lc rgb darkc(m) notitle
}
