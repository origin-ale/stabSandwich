set term png size 850, 600
set output "output/fSimDynamics_graph.png"

set key box width 2 height 1
set title "T-doped fSim magnetization dynamics" font ",24"

set label "N=10, p=.5" at graph .95,.75 right font ",20"

set xlabel "Gate" font ",16"
set ylabel "Magnetization" font ",16"
set key font ",16"
plot "output/fSimDynamics.txt" index 0 with linespoints title "CAMPS-PP (χ = 64, Nmax = 200)" lc rgb "blue", \
     "output/fSimDynamics.txt" index 1 with linespoints title "CAMPS (χ = 64)" lc rgb "red", \
     "output/fSimDynamics.txt" index 2 with linespoints title "Pauli propagation (Nmax = 200)" lc rgb "dark-green"