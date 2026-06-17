N="12"
Nsamples="100"
magic_prob="1.0"
magic_phase="pi8"

set term png size 850, 600
outfile = "output/TMD_".N."_".Nsamples."_graph.png"
set output outfile

set key bottom

set title "XXZ Floquet dynamics" font ",24"

set label "N = ".N.", Δ = 1, ".magic_phase." doping\n".Nsamples." samples" at graph .05,.10 left font ",20"

set logscale xy

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"

datafile = "output/TMD_".N."_".Nsamples.".txt"

# Whitespace-separated list of block indices to plot (0-based: index 0 is the
# first data block).
blocks = "1 3 5 7"

# Power-law fit a*x**b for each block: 1 = overlay the fit lines, 0 = only
# print the fitted parameters to the terminal.
plot_fit = 0

# Title for block `index i`: pulls the matching "# μ = …, p = …" header from the
# file (data index i corresponds to the (i+1)-th header line).
blocktitle(i) = system(sprintf("grep -a '^# .*p =' '%s' | sed -n '%dp' | sed 's/^# //'", datafile, i + 1))

# Power-law fits a*x**b for each dataset, without printing the fit log
set fit quiet
set fit logfile "/dev/null"
set fit errorvariables

# `fit` can't run inside a `plot for`, so fit every block first and stash the
# parameters in arrays keyed by the loop index.
nb = words(blocks)
array A[nb]
array B[nb]
do for [j=1:nb] {
  i = int(word(blocks, j))
  a = 1; b = 0.66
  fit [1:*] a*x**b datafile index i using 1:2 via a, b
  A[j] = a; B[j] = b
  print sprintf("%s:  a = %.4g +/- %.2g,  b = %.4g +/- %.2g", \
                blocktitle(i), a, a_err, b, b_err)
}

if (plot_fit) {
  plot for [j=1:nb] datafile index int(word(blocks, j)) \
         with yerrorlines lc j title blocktitle(int(word(blocks, j))), \
       for [j=1:nb] A[j]*x**B[j] lc j dt 2 notitle, \
       2*x lc rgb "gray" title "2t", \
       x**0.66 lc rgb "light-gray" notitle "t^{2/3}"
} else {
  plot for [j=1:nb] datafile index int(word(blocks, j)) \
         with yerrorlines lc j title blocktitle(int(word(blocks, j))), \
       2*x lc rgb "gray" title "2t", \
       x**0.66 lc rgb "light-gray" notitle "t^{2/3}"
}

print "Plotted ".outfile