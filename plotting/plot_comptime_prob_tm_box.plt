# Box plots of the per-sample computation time against magic probability.
#
# Input: output/*/comptime_prob_tm.txt (and gctime_prob_tm.txt) as written by
# Comptime_prob_tm.jl, with seven columns per method:
#   1 p, then per method  mean  std.err  q1  median  q3  min  max
# so method j (0-based, in the `methods` order of the run) starts at column
# 2 + 7*j. The box spans the quartiles, the whiskers the extrema over samples.
#
# This is the distribution behind plot_comptime_prob_tm.plt, which shows only
# the mean and its standard error.

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

campspp = 'output/carlos_varseed/comptime_prob_tm.txt'
pp = 'output/carlos_varseed_pp/comptime_prob_tm.txt'

# == Wall-clock time ====================================================================
set output "output/CompTimes_prob_tm_box.png"

plot campspp using ($1-dx):col(0,Q1):col(0,MIN):col(0,MAX):col(0,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#4269d0" \
       title "CAMPS-PP, thl 10^{-10}", \
     campspp using ($1-dx):col(0,MED):col(0,MED):col(0,MED):col(0,MED) \
       with candlesticks lw 3 lc rgb "#1b3a8c" notitle, \
     campspp using ($1-dx):col(0,MEAN) with points pt 7 ps 0.8 \
       lc rgb "#1b3a8c" title "CAMPS-PP mean", \
     pp using ($1+dx):col(0,Q1):col(0,MIN):col(0,MAX):col(0,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#efb118" \
       title "PP, thl 10^{-10}", \
     pp using ($1+dx):col(0,MED):col(0,MED):col(0,MED):col(0,MED) \
       with candlesticks lw 3 lc rgb "#9c7016" notitle, \
     pp using ($1+dx):col(0,MEAN) with points pt 7 ps 0.8 \
       lc rgb "#9c7016" title "PP mean"

# == GC time ============================================================================
set output "output/GCTimes_prob_tm_box.png"
set ylabel "GC time (s)"

campspp_gc = 'output/carlos_varseed/gctime_prob_tm.txt'
pp_gc = 'output/carlos_varseed_pp/gctime_prob_tm.txt'

plot campspp_gc using ($1-dx):col(0,Q1):col(0,MIN):col(0,MAX):col(0,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#4269d0" \
       title "CAMPS-PP, thl 10^{-10}", \
     campspp_gc using ($1-dx):col(0,MED):col(0,MED):col(0,MED):col(0,MED) \
       with candlesticks lw 3 lc rgb "#1b3a8c" notitle, \
     pp_gc using ($1+dx):col(0,Q1):col(0,MIN):col(0,MAX):col(0,Q3) \
       with candlesticks whiskerbars 0.5 lc rgb "#efb118" \
       title "PP, thl 10^{-10}", \
     pp_gc using ($1+dx):col(0,MED):col(0,MED):col(0,MED):col(0,MED) \
       with candlesticks lw 3 lc rgb "#9c7016" notitle
