N="10"
Nsamples="30"
magic_prob="1"

set term png size 850, 600
set output "output/TMD_".N."_".magic_prob."_".Nsamples."_graph.png"

set key box width 2 height 1 left

set title "p(π/8)=".magic_prob." XXZ Floquet dynamics" font ",24"

set label "N = ".N." Δ = 1, ".Nsamples." samples" at graph .05,.65 left font ",20"

set logscale xy

array mus[4] = ["0.3", "0.6", "1", "10"]

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"
plot for [i = 0:3] "output/TMD_".N."_".magic_prob."_".Nsamples.".txt" index i with yerrorlines title "µ = ".mus[i+1] ,\
x**0.66 title "t^{2/3}" lc rgb "dark-gray",\
2*x notitle lc rgb "red"
print "Plotted output/TMD_".N."_".magic_prob."_".Nsamples.".txt"
