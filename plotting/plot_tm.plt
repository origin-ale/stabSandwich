set term png size 850, 600
set output "output/TMDynamics_graph.png"

set key box width 2 height 1

set title "XXZ Floquet dynamics" font ",24"

# set label "N = 22, p=0, Δ = 1, µ = 0.8" at graph .05,.75 left font ",20"

set logscale xy
set xrange [1:12]
a = 1.0
gamma = 0.66
f(x) = a*x**gamma
fit f(x) "output/TMDynamics.txt" index 0 using 1:2 via a, gamma
gamma_str = sprintf("%.2f", gamma)

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"
plot "output/TMDynamics.txt" index 0 with yerrorlines title "CAMPS-PP" lc rgb "red", \
f(x) title "Power law (γ = ".gamma_str." )" lc rgb "purple", \
"output/TMDynamics.txt" index 1 with yerrorlines title "Pauli propagation" lc rgb "dark-green"
# 2*x title "Upper bound 2t" lc rgb "gray"