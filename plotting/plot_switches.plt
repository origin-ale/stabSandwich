# CAMPS-PP computation time against the switch point c, one curve per depth M.
#
# Input: M<M>_switch_optimization_avgs_<N>.txt as written by switch_optimization.jl.
# dir  input directory (default: output/switch_optimization, where the script writes)
# N    system size (default "20")
# Ms   the M values to plot, space-separated (default "6 9 12 15 18")
# Old data: gnuplot -e "dir='output/old/carlos_switch'" plotting/plot_switches.plt
if (!exists("dir")) { dir = "output/switch_optimization" }
if (!exists("N")) { N = "20" }
if (!exists("Ms")) { Ms = "6 9 12 15 18" }

set term png size 850, 600 font ",16"

set output "output/varM_switch_optimization_graph_total.png"

set title "CAMPS-PP computation times, ".N."×(".N."+M) circuit" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

set xrange [*:14.5]
set yrange [*:30]
set key bottom

plot for [M in Ms] dir.'/M'.M.'_switch_optimization_avgs_'.N.'.txt' using 1:2:3 with yerrorlines title "M=".M
# Old deep runs (N = 50 and 100, switch at t = 62-c and 112-c):
#   output/old/carlos_switch/deep_switch_optimization_avgs_{50,100}.txt
