# CAMPS-PP computation time against the switch point c for the N×N circuit
# (M = 0), one curve per N.
#
# Input: M0_switch_optimization_avgs_<N>.txt as written by switch_optimization.jl
# with M=0 and an explicit switch-point list, e.g.
#   julia --project=. switch_optimization.jl N=20 M=0 switchpoints=0,1,2,3,4,5,6,7,8,9,10,11,12 outdir=output/switch_optimization_square
# (a folder of their own, so that extract_minima.sh does not mix them into the M scan)
# dir  input directory (default: output/switch_optimization_square)
# Ns   the N values to plot, space-separated (default "20 50 80 100")
# Old data: gnuplot -e "dir='output/old/carlos_switch_square'" plotting/plot_switches_square.plt
if (!exists("dir")) { dir = "output/switch_optimization_square" }
if (!exists("Ns")) { Ns = "20 50 80 100" }

set term png size 850, 600 font ",16"

set output "output/square_switch_optimization_graph_total.png"

set title "CAMPS-PP computation times, NxN circuit, switch at N-c" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

set xrange [*:14.5]
set yrange [*:30]
set key bottom

plot for [N in Ns] dir.'/M0_switch_optimization_avgs_'.N.'.txt' using 1:2:3 with yerrorlines title "N=".N
