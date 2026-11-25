# CAMPS-PP computation time against the switch point c, one curve per depth M (or per N, see runs).
#
# Input: M<M>_switch_optimization_avgs_<N>.txt as written by switch_optimization.jl.
# dir  input directory (default: the newest fetched run's copy of output/switch_optimization;
#      run=<id> picks another, see run.plt). Plots go to output/plots/<run id>.
# runs several fetched runs, space-separated run ids, one curve per avgs file in each,
#      any N and M (N and Ms then filter; axes autoscale). Curves are labelled by
#      whichever of N, M varies. Plots go to output/plots/<id1>+<id2>+…
# N    system size (default "20"; with runs, all)
# Ms   the M values to plot, space-separated (default "6 9 12 15 18"; with runs, all)
# Old data: gnuplot -e "dir='output/old/carlos_switch'" plotting/plot_switches.plt
# Several runs: gnuplot -e "runs='20261123-091438-switch_optimization 20261124-141319-switch_optimization'" plotting/plot_switches.plt
datadir = "output/switch_optimization"
load "plotting/run.plt"
avgsfile(d, M, n) = d.'/M'.M.'_switch_optimization_avgs_'.n.'.txt'
# files: the avgs files to plot; Mof, Nof: their M and N values, in the same order
files = ""
Mof = ""
Nof = ""
if (exists("runs")) {
  dirs = ""
  do for [r in runs] {
    d = system(sprintf("plotting/run_dirs.sh data '%s' '%s'", datadir, r))
    if (strlen(d) == 0) { exit status 1 }
    dirs = dirs." ".d
    # "M,N" pairs, by N then M
    found = system(sprintf("ls '%s' | sed -n 's/^M\\([0-9]*\\)_switch_optimization_avgs_\\([0-9]*\\)\\.txt$/\\1,\\2/p' | sort -t, -k2n -k1n | paste -sd ' ' -", d))
    do for [mn in found] {
      M = mn[1:strstrt(mn, ",") - 1]
      n = mn[strstrt(mn, ",") + 1:]
      if ((!exists("Ms") || strstrt(" ".Ms." ", " ".M." ") > 0) && (!exists("N") || n eq N)) {
        files = files." ".avgsfile(d, M, n); Mof = Mof." ".M; Nof = Nof." ".n
      }
    }
  }
  plotdir = plot_dir(dirs)
} else {
  if (!exists("dir")) { dir = rundir }
  plotdir = plot_dir(dir)
  if (!exists("N")) { N = "20" }
  if (!exists("Ms")) { Ms = "6 9 12 15 18" }
  do for [M in Ms] { files = files." ".avgsfile(dir, M, N); Mof = Mof." ".M; Nof = Nof." ".N }
}
if (words(files) == 0) {
  print "plot_switches.plt: no M*_switch_optimization_avgs_*.txt to plot"
  exit status 1
}
ndistinct(l) = system("echo '".l."' | tr ' ' '\\n' | grep . | sort -u | wc -l") + 0
varN = ndistinct(Nof) > 1
varM = ndistinct(Mof) > 1 || !varN
N1 = word(Nof, 1)
M1 = word(Mof, 1)
label(k) = (varN ? "N=".word(Nof, k) : "").(varN && varM ? ", " : "").(varM ? "M=".word(Mof, k) : "")

set term png size 850, 600 font ",16"

set output plotdir."/".(varN ? "varN" : "").(varM ? "varM" : "")."_switch_optimization_graph_total".(varN && !varM ? "_M".M1 : "").".png"

set title "CAMPS-PP computation times, ".(varN ? "N×(N+".(varM ? "M" : M1).")" : N1."×(".N1."+M)")." circuit" font ",24"
set logscale y
set format y "%.1f"
set xlabel "c"
set ylabel "Time (s)"

if (exists("runs")) {
  set format y "%1g"
  set xtics 1
  set grid ytics
  set offsets 0.5, 0.5, 0, 0
  set key top right
} else {
  set xrange [*:14.5]
  set yrange [*:30]
  set key bottom
}

plot for [k = 1:words(files)] word(files, k) using 1:2:3 with yerrorlines title label(k)
# Old deep runs (N = 50 and 100, switch at t = 62-c and 112-c):
#   output/old/carlos_switch/deep_switch_optimization_avgs_{50,100}.txt
