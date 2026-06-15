N="12"
Nsamples="100"
magic_prob="1.0"

set term png size 850, 600
set output "output/TMD_".N."_".magic_prob."_".Nsamples."_graph.png"

set key bottom

set title "Dope=".magic_prob." XXZ Floquet dynamics" font ",24"

set label "N = ".N.", Δ = 1\n".Nsamples." samples" at graph .05,.10 left font ",20"

set logscale xy

array mus[4] = ["0.3", "0.6", "1", "10"]

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"

datafile = "output/TMD_".N."_".magic_prob."_".Nsamples.".txt"

# Power-law fits a*x**b for each dataset, without printing the fit log
set fit quiet
set fit logfile "/dev/null"
set fit errorvariables

array A[4]
array B[4]
array Aerr[4]
array Berr[4]
do for [i = 0:3] {
    a = 1; b = 0.66
    fit [1:*] a*x**b datafile index i using 1:2 via a, b
    A[i+1] = a; B[i+1] = b; Aerr[i+1] = a_err; Berr[i+1] = b_err
}
print "Power-law fits  f(x) = a*x**b"
do for [i = 0:3] {
    print sprintf("  µ = %-4s  a = %.4f +/- %.4f   b = %.4f +/- %.4f", \
        mus[i+1], A[i+1], Aerr[i+1], B[i+1], Berr[i+1])
}

plot for [i = 0:3] datafile index i with yerrorlines title "µ = ".mus[i+1] ,\
x**0.66 title "t^{2/3}" lc rgb "dark-gray",\
2*x title "2t (upper bound)" lc rgb "red"
print "Plotted output/TMD_".N."_".magic_prob."_".Nsamples.".txt"
