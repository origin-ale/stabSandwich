set term png size 850, 600
set output "output/minima_graph.png"

set title "Optimal-switch computation time, N=20" font ",24"

set logscale y
set xrange [5.5:18.5]

set xlabel "M" font ",16"
set ylabel "Time (s)" font ",16"
set key font ",16"

f(x) = c + a*2**(b*x**2)
a = 0.1
b = 0.015
c = 0.2
fit f(x) "output/switch/switch_minima.txt" using 1:3:4 yerrors via a,b,c

plot "output/switch/switch_minima.txt" using 1:3:4 with yerrorlines notitle,\
f(x) notitle lc rgb "gray"