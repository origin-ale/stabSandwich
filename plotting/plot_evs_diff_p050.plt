set term png size 850, 600 font ",16"

set xlabel "Layer"
set xrange [0:*]

input1 = 'output/carlos_varseed/evs_prob_tm_campspp_10.txt'
input2 = 'output/carlos_varseed_pp/evs_prob_tm_pp_10.txt'

# p blocks are separated by double blank lines: p = 0.055 is index 11
idx = 10

# paste joins the two files line by line (identical block structure),
# so file 2's columns become 4:5:6
paired = '< paste '.input1.' '.input2

# --- expectation values ---
set output "output/EVs_p050.png"
set title "46-qubit 23-layer XXZ circuit, p = 0.050" font ",24"
set ylabel "Expectation value"
set key top left

plot input1 index idx using 1:2:3 with yerrorlines title "CAMPS-PP (thl 10^{-10})",\
     input2 index idx using 1:2:3 with yerrorlines title "PP (thl 10^{-10})"

# --- pointwise difference ---
set output "output/EVs_p050_diff.png"
set title "46-qubit 23-layer XXZ circuit, p = 0.050" font ",24"
set ylabel "Difference"
set key top left

plot paired index idx using 1:($2-$5) with linespoints title "diff CAMPS-PP - PP, thl 10^{-10}"
