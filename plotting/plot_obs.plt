set term png
set output "output/CircuitDynamics_graph.png"

set key box width 2 height 1
set title "Magnetization dynamics, N = 13"
set xlabel "Time"
set ylabel "Magnetization"
plot "output/CircuitDynamics.txt" index 0 with linespoints title "CAMPS-PP" lc rgb "blue", \
     "output/CircuitDynamics.txt" index 1 with linespoints title "CAMPS" lc rgb "red", \
     "output/CircuitDynamics.txt" index 2 with linespoints title "Pauli Prop." lc rgb "dark-green"