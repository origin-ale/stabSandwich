set term png size 850, 600 font ",16"
set output "output/EVs_p055_carlos10_vs_carlos9.png"

set title "46-qubit 23-layer XXZ circuit, p = 0.055" font ",24"
set xlabel "Layer"
set ylabel "Expectation value"
set y2label "Difference"

set ytics nomirror
set y2tics

set xrange [0:*]

set key top left

input10 = 'output/carlos10/evs_prob_tm_campspp_10.txt'
input9  = 'output/carlos5/evs_prob_tm_campspp_5.txt'

# p blocks are separated by double blank lines: p = 0.055 is index 11
idx = 11

# paste joins the two files line by line (identical block structure),
# so file 9's columns become 4:5:6
paired = '< paste '.input10.' '.input9

plot input10 index idx using 1:2:3 with yerrorlines title "CAMPS-PP (thl 10^{-10})",\
     input9  index idx using 1:2:3 with yerrorlines title "CAMPS-PP (thl 10^{-5})",\
     paired  index idx using 1:($2-$5) axes x1y2 with linespoints title "difference"
