# Box plots of the per-sample computation time against magic probability.
#
# Input: comptime_prob_tm.txt (and gctime_prob_tm.txt) as written by
# Comptime_prob_tm.jl, with seven columns per method:
#   1 p, then per method  mean  std.err  q1  median  q3  min  max
# so method j (0-based, in the `methods` order of the run) starts at column
# 2 + 7*j. The box spans the quartiles, the whiskers the extrema over samples.
#
# This is the distribution behind plot_comptime_prob_tm.plt, which shows only
# the mean and its standard error.
#
# dir_campspp, dir_pp  directories holding each method's files
#                      (default: output/Comptime_prob_tm for both, i.e. one run
#                      with methods=campspp,pp)
# j_campspp, j_pp      each method's 0-based position in its file's `methods`
#                      (default 0 and 1)
# Old data (one method per run):
#   gnuplot -e "dir_campspp='output/old/carlos_indep'; dir_pp='output/old/carlos_indep_pp'; j_pp=0" plotting/plot_comptime_prob_tm_box.plt
if (!exists("dir_campspp")) { dir_campspp = "output/Comptime_prob_tm" }
if (!exists("dir_pp")) { dir_pp = "output/Comptime_prob_tm" }
if (!exists("j_campspp")) { j_campspp = 0 }
if (!exists("j_pp")) { j_pp = 1 }

set term png size 850, 600 font ",16"

set title "46-qubit 23-layer XXZ circuit" font ",24"
set logscale y
set format y "%1g"
set xlabel "Magic probability"
set ylabel "Time (s)"

set xrange [*:*]
set yrange [*:*]
set key top left

set style fill solid 0.25 border lc rgb "black"

# Column of quantity `q` for method j (0-based) in a file with 7 columns/method.
col(j, q) = 2 + 7*j + q
Q1 = 2; MED = 3; Q3 = 4; MIN = 5; MAX = 6; MEAN = 0

# The p spacing is 0.005, so boxes must be far narrower than that; `dx` also
# offsets the two methods sideways so their boxes sit side by side instead of
# on top of each other.
bw = 0.0012
set boxwidth bw absolute
dx = 0.0009

campspp = dir_campspp.'/comptime_prob_tm.txt'
pp = dir_pp.'/comptime_prob_tm.txt'

# == Wall-clock time ====================================================================
set output "output/CompTimes_prob_tm_box.png"

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
set output "output/GCTimes_prob_tm_box.png"
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
