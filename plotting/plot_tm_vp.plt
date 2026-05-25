set term png size 850, 600
set output "output/TMDynamics_vp_graph.png"

set key box width 2 height 1 left

set title "π/8 XXZ Floquet dynamics" font ",24"

set label "N = 12, Δ = 1, 25 samples" at graph .05,.65 left font ",20"

set logscale xy
set xrange [.9:6]
set linetype 5 lc rgb "dark-gray"

array mus[4] = ["0.3", "0.6", "1", "10"]

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"
plot for [i = 0:3] "output/TMDynamics_varparams.txt" index i with yerrorlines title "µ = ".mus[i+1] ,\
x**0.66 title "t^{2/3}"