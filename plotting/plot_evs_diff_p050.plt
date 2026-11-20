# Mean expectation value per layer in one p block, CAMPS-PP against PP, and their
# pointwise difference.
#
# Input: evs_prob_tm_{campspp,pp}_<thl exponent>.txt as written by Comptime_prob_tm.jl.
# dir_campspp, dir_pp  directories holding each method's file
#                      (default: output/Comptime_prob_tm for both)
# thl_campspp, thl_pp  exponent suffix of each file name (default "10", i.e. thl 1e-10)
# idx                  0-based p block to plot (default 10, which is p = 0.050 in the
#                      old runs: p from 0 in steps of 0.005)
# Old data:
#   gnuplot -e "dir_campspp='output/old/carlos_varseed'; dir_pp='output/old/carlos_varseed_pp'" plotting/plot_evs_diff_p050.plt
if (!exists("dir_campspp")) { dir_campspp = "output/Comptime_prob_tm" }
if (!exists("dir_pp")) { dir_pp = "output/Comptime_prob_tm" }
if (!exists("thl_campspp")) { thl_campspp = "10" }
if (!exists("thl_pp")) { thl_pp = "10" }
if (!exists("idx")) { idx = 10 }

set term png size 850, 600 font ",16"

set xlabel "Layer"
set xrange [0:*]

input1 = dir_campspp.'/evs_prob_tm_campspp_'.thl_campspp.'.txt'
input2 = dir_pp.'/evs_prob_tm_pp_'.thl_pp.'.txt'

# "p = …" of block idx, from the file's "# n = …, μ = …" header ("# p = …" in p-based data) lines
ptitle = system(sprintf("grep -a '^# [pn] =' '%s' | sed -n '%dp' | sed 's/^# //; s/,.*//'", input1, idx + 1))

# paste joins the two files line by line (identical block structure),
# so file 2's columns become 4:5:6
paired = '< paste '.input1.' '.input2

# --- expectation values ---
set output "output/EVs_p050.png"
set title "46-qubit 23-layer XXZ circuit, ".ptitle font ",24"
set ylabel "Expectation value"
set key top left

plot input1 index idx using 1:2:3 with yerrorlines title "CAMPS-PP (thl 10^{-".thl_campspp."})",\
     input2 index idx using 1:2:3 with yerrorlines title "PP (thl 10^{-".thl_pp."})"

# --- pointwise difference ---
set output "output/EVs_p050_diff.png"
set title "46-qubit 23-layer XXZ circuit, ".ptitle font ",24"
set ylabel "Difference"
set key top left

plot paired index idx using 1:($2-$5) with linespoints title "diff CAMPS-PP - PP, thl 10^{-".thl_campspp."}"
