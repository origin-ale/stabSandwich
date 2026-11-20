# Transferred magnetization against time, one image per n value (all μ of that n;
# p in p-based data).
#
# Input: <prefix>_<N>_<Nsamples>.txt as written by TransferredMag_varparams.jl.
# dir         input directory (default: output/TransferredMag_varparams, where the
#             script writes)
# N, Nsamples, prefix  as in the file name (default "46", "50", "TMDxy")
# plot_fit    1 = overlay the power-law fits, 0 = only print them (default 0)
# Old data: gnuplot -e "dir='output/old/TMD'; N='12'; Nsamples='100'" plotting/plot_tm_vp_v2.plt
# (the 12- and 46-qubit runs are also in output/carlos_tmd{12,46}{xy,z}).
if (!exists("dir")) { dir = "output/TransferredMag_varparams" }
if (!exists("N")) { N = "46" }
if (!exists("Nsamples")) { Nsamples = "50" }
if (!exists("prefix")) { prefix = "TMDxy" }
if (!exists("plot_fit")) { plot_fit = 0 }

set term png size 850, 600

set key bottom

set title "XXZ Floquet dynamics" font ",24"

set label "N = ".N.", ϕ = θ = π/4\nDoping with 3π/16 on XX,YY\n".Nsamples." samples" at graph .05,.15 left font ",16"

set logscale xy
set yrange [0.01:*]

set xlabel "Time" font ",16"
set ylabel "Transferred magnetization" font ",16"
set key font ",16"

datafile = dir."/".prefix."_".N."_".Nsamples.".txt"

# Blocks are written μ-major ("# n = …, μ = …", n varying fastest; p in p-based data), so with np
# distinct n values, image k shows blocks k-1, k-1+np, k-1+2np, … (one n, every μ).
# Each entry of `blocksets` is that whitespace-separated list of 0-based indices.
nb_all = system(sprintf("grep -ac '^# [pn] =' '%s'", datafile)) + 0
np = system(sprintf("grep -a '^# [pn] =' '%s' | sed 's/^# [pn] = //; s/,.*//' | sort -u | wc -l", datafile)) + 0
array blocksets[np]
do for [k=1:np] {
  blocksets[k] = ""
  do for [i=k-1:nb_all-1:np] { blocksets[k] = blocksets[k].sprintf(" %d", i) }
}

# Fit window: points with time >= fitmin; a block needs at least 3 there
# (2 parameters), or its fit is skipped.
fitmin = 10
nfit(i) = system(sprintf("awk 'BEGIN{b=-1} /^# [pn] =/{b++; next} /^#/{next} b==%d && NF>=2 && $1>=%g {c++} END{print c+0}' '%s'", i, fitmin, datafile)) + 0

# Title for block `index i`: pulls the matching "# n = …, μ = …" header ("# p = …" in p-based data) from the
# file (data index i corresponds to the (i+1)-th header line).
blocktitle(i) = system(sprintf("grep -a '^# [pn] =' '%s' | sed -n '%dp' | sed 's/^# //'", datafile, i + 1))

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
    if (nfit(i) >= 3) {
      a = 1; b = 0.66
      fit [fitmin:*] a*x**b datafile index i using 1:2 via a, b
      A[j] = a; B[j] = b
      print sprintf("%s:  a = %.4g +/- %.2g,  b = %.4g +/- %.2g", \
                    blocktitle(i), a, a_err, b, b_err)
    } else {
      A[j] = NaN; B[j] = NaN
      print sprintf("%s:  fit skipped (fewer than 3 points at time >= %g)", blocktitle(i), fitmin)
    }
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