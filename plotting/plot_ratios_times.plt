# Mean computation times against t/N, exact CAMPS against exact PP.
#
# Input: comptimes_ratios_<N>_avgs.txt as written by Comptime_ratios.jl.
# dir  input directory (default: the newest fetched run's copy of output/Comptime_ratios;
#      run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# N    system size in the file name (default "40")
# Old data: gnuplot -e "dir='output/old/carlos_ratios'" plotting/plot_ratios_times.plt
datadir = "output/Comptime_ratios"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)
if (!exists("N")) { N = "40" }

set term png size 850, 600 font ",16"

set output plotdir."/CompTimes_ratios_".N."_graph.png"

set title "Computation times, N=".N." t-rotation circuit" font ",24"
set xrange [0.65:1.25]
set logscale y
set format y "%.1e"
set xlabel "t/N"
set ylabel "Time (s)"

set key bottom

avgs = dir.'/comptimes_ratios_'.N.'_avgs.txt'
plot avgs using 1:2:3 with yerrorlines title "Exact CAMPS",\
     avgs using 1:4:5 with yerrorlines title "Exact PP"
