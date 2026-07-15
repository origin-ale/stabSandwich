set term png size 850, 600 font ",16"
set output "output/Resources_prob_tm_bd_camps_graph.png"

datafile = "output/resources_prob_tm_bd_camps.txt"

set title "46-qubit XXZ circuit — CAMPS bond dimension" font ",24"
set xlabel "Gate" font ",16"
set ylabel "Mean bond dimension" font ",16"
set key font ",16" top left

# Title for data index i: pulls the matching "# p = …, μ = …" header line
# from the file (data index i corresponds to the (i+1)-th such header line).
blocktitle(i) = system(sprintf("grep -a '^# p =' '%s' | sed -n '%dp' | sed 's/^# //'", datafile, i + 1))

nb = system(sprintf("grep -ac '^# p =' '%s'", datafile)) + 0

plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title blocktitle(i)
