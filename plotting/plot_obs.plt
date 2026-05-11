set term png
set output "output/NaiveDynamics_graph.png"

set key box width 2 height 1
set xlabel "Time"
set ylabel "[Observable]"
plot "output/NaiveDynamics.txt" with linespoints title "CAMPS-PP"