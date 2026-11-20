# Mean expectation value per layer in one n block (p block in p-based data),
# CAMPS-PP against PP, and their pointwise difference.
#
# Input: evs_prob_tm_{campspp,pp}_<thl exponent>.txt as written by Comptime_prob_tm.jl.
# dir_campspp, dir_pp  directories holding each method's file
#                      (default for both: the newest fetched run's copy of
#                      output/Comptime_prob_tm; run=<id> picks another, see run.plt)
# thl_campspp, thl_pp  exponent suffix of each file name (default "10", i.e. thl 1e-10)
# idx                  0-based block to plot (default: the last block both files
#                      have, i.e. the most doped n or the largest p)
# Output: EVs_<block>.png and EVs_<block>_diff.png in output/plots/<run id>, where <block> is
# the block's label, e.g. EVs_n45.png or EVs_p0.05.png.
# Old data (p = 0.05 is block 10: p from 0 in steps of 0.005):
#   gnuplot -e "dir_campspp='output/old/carlos_varseed'; dir_pp='output/old/carlos_varseed_pp'; idx=10" plotting/plot_evs_diff_p050.plt
datadir = "output/Comptime_prob_tm"
load "plotting/run.plt"
if (!exists("dir_campspp")) { dir_campspp = rundir }
if (!exists("dir_pp")) { dir_pp = rundir }
plotdir = plot_dir(dir_campspp." ".dir_pp)
if (!exists("thl_campspp")) { thl_campspp = "10" }
if (!exists("thl_pp")) { thl_pp = "10" }

set term png size 850, 600 font ",16"

set xlabel "Layer"
set xrange [0:*]

input1 = dir_campspp.'/evs_prob_tm_campspp_'.thl_campspp.'.txt'
input2 = dir_pp.'/evs_prob_tm_pp_'.thl_pp.'.txt'

# Blocks present in both files (the old runs have different counts per method)
nblocks(f) = system(sprintf("grep -ac '^# [pn] =' '%s'", f)) + 0
nb = (nblocks(input1) < nblocks(input2)) ? nblocks(input1) : nblocks(input2)
if (!exists("idx")) { idx = nb - 1 }
if (idx < 0 || idx >= nb) {
  print sprintf("plot_evs_diff_p050.plt: idx=%d out of range, the files share blocks 0..%d", idx, nb - 1)
  exit status 1
}

# "n = …" (or "p = …") of block idx, from the file's "# n = …, μ = …" header lines
ptitle = system(sprintf("grep -a '^# [pn] =' '%s' | sed -n '%dp' | sed 's/^# //; s/,.*//'", input1, idx + 1))
# File-name tag: "n = 45" -> "n45", "p = 0.05" -> "p0.05"
tag = system(sprintf("echo '%s' | tr -d ' ='", ptitle))

# paste joins the two files line by line (identical block structure),
# so file 2's columns become 4:5:6
paired = '< paste '.input1.' '.input2

# --- expectation values ---
set output plotdir."/EVs_".tag.".png"
set title "46-qubit 23-layer XXZ circuit, ".ptitle font ",24"
set ylabel "Expectation value"
set key top left

plot input1 index idx using 1:2:3 with yerrorlines title "CAMPS-PP (thl 10^{-".thl_campspp."})",\
     input2 index idx using 1:2:3 with yerrorlines title "PP (thl 10^{-".thl_pp."})"

# --- pointwise difference ---
set output plotdir."/EVs_".tag."_diff.png"
set title "46-qubit 23-layer XXZ circuit, ".ptitle font ",24"
set ylabel "Difference"
set key top left

plot paired index idx using 1:($2-$5) with linespoints title "diff CAMPS-PP - PP, thl 10^{-".thl_campspp."}"
