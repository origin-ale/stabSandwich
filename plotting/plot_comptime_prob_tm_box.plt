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
# dir_campspp, dir_pp  directories holding each method's files
#                      (default for both: the newest fetched run's copy of
#                      output/Comptime_prob_tm, i.e. one run with methods=campspp,pp;
#                      run=<id> picks another, see run.plt)
# Plots go to output/plots/<run id> (<id1>+<id2> for dirs from two runs).
# j_campspp, j_pp      each method's 0-based position in its file's `methods`
#                      (default 0 and 1)
# Old data (one method per run):
#   gnuplot -e "dir_campspp='output/old/carlos_indep'; dir_pp='output/old/carlos_indep_pp'; j_pp=0" plotting/plot_comptime_prob_tm_box.plt
datadir = "output/Comptime_prob_tm"
load "plotting/run.plt"
if (!exists("dir_campspp")) { dir_campspp = rundir }
if (!exists("dir_pp")) { dir_pp = rundir }
plotdir = plot_dir(dir_campspp." ".dir_pp)
if (!exists("j_campspp")) { j_campspp = 0 }
if (!exists("j_pp")) { j_pp = 1 }
# x axis: n (doped triples per sample) for count-mode runs, p for p-based data
xlab(f) = (system(sprintf("grep -ac 'magic_counts=' '%s'", f)) + 0 > 0) ? "Doped triples per sample" : "Magic probability"

set term png size 850, 600 font ",16"

set title "46-qubit 23-layer XXZ circuit" font ",24"
set logscale y
set format y "%1g"
set xlabel xlab(dir_campspp.'/comptime_prob_tm.txt')
set ylabel "Time (s)"

set xrange [*:*]
set yrange [*:*]
set key top left

set style fill solid 0.25 border lc rgb "black"

# Column of quantity `q` for method j (0-based) in a file with 7 columns/method.
col(j, q) = 2 + 7*j + q
Q1 = 2; MED = 3; Q3 = 4; MIN = 5; MAX = 6; MEAN = 0

campspp = dir_campspp.'/comptime_prob_tm.txt'
pp = dir_pp.'/comptime_prob_tm.txt'

# The p spacing is 0.005, so boxes must be far narrower than that; `dx` also
# offsets the two methods sideways so their boxes sit side by side instead of
# on top of each other. For count-mode runs both scale with the smallest gap
# between n values (the same fractions: 0.24 and 0.18 of it).
bw = 0.0012
dx = 0.0009
if (system(sprintf("grep -ac 'magic_counts=' '%s'", campspp)) + 0 > 0) {
  xstep = system(sprintf("awk '!/^#/ && NF {x[k++]=$1} END {m=0; for (i=1;i<k;i++) {d=x[i]-x[i-1]; if (d>0 && (m==0 || d<m)) m=d}; print (m>0 ? m : 1)}' '%s'", campspp)) + 0
  bw = 0.24*xstep
  dx = 0.18*xstep
}
set boxwidth bw absolute

# == Wall-clock time ====================================================================
set output plotdir."/CompTimes_prob_tm_box.png"

plot campspp using ($1-dx):col(j_campspp,Q1):col(j_campspp,MIN):col(j_campspp,MAX):col(j_campspp,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#4269d0" \
       title "CAMPS-PP, thl 10^{-10}", \
     campspp using ($1-dx):col(j_campspp,MED):col(j_campspp,MED):col(j_campspp,MED):col(j_campspp,MED) \
       with candlesticks lw 3 lc rgb "#1b3a8c" notitle, \
     campspp using ($1-dx):col(j_campspp,MEAN) with linespoints pt 7 ps 0.8 \
       lc rgb "#1b3a8c" title "CAMPS-PP mean", \
     pp using ($1+dx):col(j_pp,Q1):col(j_pp,MIN):col(j_pp,MAX):col(j_pp,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#efb118" \
       title "PP, thl 10^{-10}", \
     pp using ($1+dx):col(j_pp,MED):col(j_pp,MED):col(j_pp,MED):col(j_pp,MED) \
       with candlesticks lw 3 lc rgb "#9c7016" notitle, \
     pp using ($1+dx):col(j_pp,MEAN) with linespoints pt 7 ps 0.8 \
       lc rgb "#9c7016" title "PP mean"

# == GC time ============================================================================
set output plotdir."/GCTimes_prob_tm_box.png"
set ylabel "GC time (s)"

campspp_gc = dir_campspp.'/gctime_prob_tm.txt'
pp_gc = dir_pp.'/gctime_prob_tm.txt'

plot campspp_gc using ($1-dx):col(j_campspp,Q1):col(j_campspp,MIN):col(j_campspp,MAX):col(j_campspp,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#4269d0" \
       title "CAMPS-PP, thl 10^{-10}", \
     campspp_gc using ($1-dx):col(j_campspp,MED):col(j_campspp,MED):col(j_campspp,MED):col(j_campspp,MED) \
       with candlesticks lw 3 lc rgb "#1b3a8c" notitle, \
     pp_gc using ($1+dx):col(j_pp,Q1):col(j_pp,MIN):col(j_pp,MAX):col(j_pp,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#efb118" \
       title "PP, thl 10^{-10}", \
     pp_gc using ($1+dx):col(j_pp,MED):col(j_pp,MED):col(j_pp,MED):col(j_pp,MED) \
       with candlesticks lw 3 lc rgb "#9c7016" notitle
