set term png size 850, 600
set output "output/CircuitDynamics_graph.png"

set key box width 2 height 1
set title "Random rotation magnetization dynamics, N = 21" font ",24"
set xlabel "Time" font ",16"
set ylabel "Magnetization" font ",16"
set key font ",16"
plot "output/CircuitDynamics.txt" index 0 with linespoints title "CAMPS-PP (χ = 64, Nmax = 200)" lc rgb "blue", \
     "output/CircuitDynamics.txt" index 1 with linespoints title "CAMPS (χ = 128)" lc rgb "red", \
     "output/CircuitDynamics.txt" index 2 with linespoints title "Pauli propagation (Nmax = 2e6)" lc rgb "dark-green"