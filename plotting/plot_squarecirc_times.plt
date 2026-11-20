# Mean computation times against N for the N-qubit N-rotation circuit.
#
# Input: comptimes_3_squarecirc_avgs.txt as written by Comptime_3_squarecirc.jl.
# dir  input directory (default: the newest fetched run's copy of output/Comptime_3_squarecirc;
#      run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# Old data: gnuplot -e "dir='output/old/carlos_squarecirc3_v2'" plotting/plot_squarecirc_times.plt
datadir = "output/Comptime_3_squarecirc"
load "plotting/run.plt"
if (!exists("dir")) { dir = rundir }
plotdir = plot_dir(dir)

set term png size 850, 600 font ",16"
set output plotdir."/CompTimes_squarecirc_graph.png"

set title "Computation times, N-qubit N-rotation circuit" font ",24"
set logscale xy
set format y "%.1e"
set xlabel "N"
set ylabel "Time (s)"

set key bottom

avgs = dir.'/comptimes_3_squarecirc_avgs.txt'
plot avgs using 1:2:3 with yerrorlines title "Exact CAMPS",\
     avgs using 1:4:5 with yerrorlines title "Exact PP",\
     avgs using 1:6:7 with yerrorlines title "Exact MPS" lc rgb "red"
