# Input: switch_minima.txt, written by plotting/extract_minima.sh into the same
# directory as the switch_optimization.jl output it reads.
# dir  input directory (default: output/switch_optimization)
# Old data: gnuplot -e "dir='output/old/carlos_switch'" plotting/plot_minima_c.plt
if (!exists("dir")) { dir = "output/switch_optimization" }

set term png size 850, 600
set output "output/minima_c_graph.png"

set title "Optimal switch position, N=20" font ",24"

# M axis: the data's M range, padded by half a step on each side
stats dir."/switch_minima.txt" using 1 nooutput
set xrange [STATS_min - 0.5 : STATS_max + 0.5]
set yrange [3:12]

set xlabel "M" font ",16"
set ylabel "c" font ",16"
set key font ",16" left top

g(x) = m*x + q
m = 0.6
q = 1.0
fit g(x) dir."/switch_minima.txt" using 1:2 via m,q

plot dir."/switch_minima.txt" using 1:2 with points pt 7 ps 1.5 title "optimal c",\
g(x) title sprintf("%.3f M + %.3f", m, q) lc rgb "gray"
