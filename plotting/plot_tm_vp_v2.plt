N="46"
Nsamples="50"
prefix="TMDxy"

set term png size 850, 600

set key bottom

set title "XXZ Floquet dynamics" font ",24"

set label "N = ".N.", ϕ = θ = π/4\nDoping with 3π/16 on XX,YY\n".Nsamples." samples" at graph .05,.15 left font ",16"

set logscale xy
set yrange [0.01:*]

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"

datafile = "output/".prefix."_".N."_".Nsamples.".txt"

# Each entry is a whitespace-separated list of block indices to plot (0-based:
# index 0 is the first data block). One image is saved per entry.
array blocksets[3]
blocksets[1] = "0 3 6 9"
blocksets[2] = "1 4 7 10"
blocksets[3] = "2 5 8 11"
# blocksets[4] = "3 7 11"

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

# Produce one image per block set.
do for [k=1:|blocksets|] {
  blocks = blocksets[k]
  outfile = "output/".prefix."_".N."_".Nsamples."_graph_".k.".png"
  set output outfile

  # `fit` can't run inside a `plot for`, so fit every block first and stash the
  # parameters in arrays keyed by the loop index.
  nb = words(blocks)
  array A[nb]
  array B[nb]
  do for [j=1:nb] {
    i = int(word(blocks, j))
    a = 1; b = 0.66
    fit [10:*] a*x**b datafile index i using 1:2 via a, b
    A[j] = a; B[j] = b
    print sprintf("%s:  a = %.4g +/- %.2g,  b = %.4g +/- %.2g", \
                  blocktitle(i), a, a_err, b, b_err)
  }

  if (plot_fit) {
    plot for [j=1:nb] datafile index int(word(blocks, j)) \
           with yerrorlines lc j title blocktitle(int(word(blocks, j))), \
         for [j=1:nb] A[j]*x**B[j] lc j dt 2 notitle, \
         2*x lc rgb "gray" title "2t", \
         x**0.66 lc rgb "light-gray" title "t^{2/3}"
  } else {
    plot for [j=1:nb] datafile index int(word(blocks, j)) \
           with yerrorlines lc j title blocktitle(int(word(blocks, j))), \
         2*x lc rgb "gray" title "2t", \
         x**0.66 lc rgb "light-gray" title "t^{2/3}"
  }

  print "Plotted ".outfile
}