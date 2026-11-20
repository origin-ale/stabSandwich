# Input: switch_minima.txt, written by plotting/extract_minima.sh into the plot
# folder output/plots/<run id> of the switch_optimization.jl output it reads.
# dir  input directory (default: the newest fetched run's copy of output/switch_optimization;
#      run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# Old data: gnuplot -e "dir='output/old/carlos_switch'" plotting/plot_minima_c.plt
datadir = "output/switch_optimization"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)

set term png size 850, 600
set output plotdir."/minima_c_graph.png"

set title "Optimal switch position, N=20" font ",24"

# M axis: the data's M range, padded by half a step on each side
stats plotdir."/switch_minima.txt" using 1 nooutput
set xrange [STATS_min - 0.5 : STATS_max + 0.5]
set yrange [3:12]

set xlabel "M" font ",16"
set ylabel "c" font ",16"
set key font ",16" left top

g(x) = m*x + q
m = 0.6
q = 1.0
fit g(x) plotdir."/switch_minima.txt" using 1:2 via m,q

plot plotdir."/switch_minima.txt" using 1:2 with points pt 7 ps 1.5 title "optimal c",\
g(x) title sprintf("%.3f M + %.3f", m, q) lc rgb "gray"
