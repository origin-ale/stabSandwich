# Input: switch_minima.txt, written by plotting/extract_minima.sh into the plot
# folder output/plots/<run id> of the switch_optimization.jl output it reads.
# dir  input directory (default: the newest fetched run's copy of output/switch_optimization;
#      run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# Old data: gnuplot -e "dir='output/old/carlos_switch'" plotting/plot_minima.plt
datadir = "output/switch_optimization"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)

set term png size 850, 600
set output plotdir."/minima_graph.png"

set title "Optimal-switch computation time, N=20" font ",24"

set logscale y
# M axis: the data's M range, padded by half a step on each side
stats plotdir."/switch_minima.txt" using 1 nooutput
set xrange [STATS_min - 0.5 : STATS_max + 0.5]

set xlabel "M" font ",16"
set ylabel "Time (s)" font ",16"
set key font ",16"

f(x) = c + a*2**(b*x**2)
a = 0.1
b = 0.015
c = 0.2
fit f(x) plotdir."/switch_minima.txt" using 1:3:4 yerrors via a,b,c

plot plotdir."/switch_minima.txt" using 1:3:4 with yerrorlines notitle,\
f(x) notitle lc rgb "gray"