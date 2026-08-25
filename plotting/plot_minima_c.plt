set term png size 850, 600
set output "output/minima_c_graph.png"

set title "Optimal switch position, N=20" font ",24"

set xrange [5.5:18.5]
set yrange [3:12]

set xlabel "M" font ",16"
set ylabel "c" font ",16"
set key font ",16" left top

g(x) = m*x + q
m = 0.6
q = 1.0
fit g(x) "output/carlos_switch/switch_minima.txt" using 1:2 via m,q

plot "output/carlos_switch/switch_minima.txt" using 1:2 with points pt 7 ps 1.5 title "optimal c",\
g(x) title sprintf("%.3f M + %.3f", m, q) lc rgb "gray"
