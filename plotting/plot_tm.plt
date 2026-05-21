set term png size 850, 600
set output "output/TMDynamics_graph.png"

set key box width 2 height 1 bottom

set title "XXZ dynamics, doping with π/8 XX" font ",24"

set label "N = 46, p=0, Δ = 1" at graph .95,.25 right font ",20"

# set logscale xy

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"
plot "output/TMDynamics.txt" index 0 with linespoints title "CAMPS-PP (χ = 64, Nmax = 200)" lc rgb "blue", \
     "output/TMDynamics.txt" index 1 with linespoints title "Pauli propagation (Nmax = 200)" lc rgb "dark-green"
    #  "output/TMDynamics.txt" index 1 with linespoints title "CAMPS (χ = 64)" lc rgb "red", \