# Input: switch_minima.txt, written by plotting/extract_minima.sh into the same
# directory as the switch_optimization.jl output it reads.
# dir  input directory (default: output/switch_optimization)
# Old data: gnuplot -e "dir='output/old/carlos_switch'" plotting/plot_minima.plt
if (!exists("dir")) { dir = "output/switch_optimization" }

set term png size 850, 600
set output "output/minima_graph.png"

set title "Optimal-switch computation time, N=20" font ",24"

set logscale y
# M axis: the data's M range, padded by half a step on each side
stats dir."/switch_minima.txt" using 1 nooutput
set xrange [STATS_min - 0.5 : STATS_max + 0.5]

set xlabel "M" font ",16"
set ylabel "Time (s)" font ",16"
set key font ",16"

f(x) = c + a*2**(b*x**2)
a = 0.1
b = 0.015
c = 0.2
fit f(x) dir."/switch_minima.txt" using 1:3:4 yerrors via a,b,c

plot dir."/switch_minima.txt" using 1:3:4 with yerrorlines notitle,\
f(x) notitle lc rgb "gray"