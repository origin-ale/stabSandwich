set term png size 850, 600 font ",16"

set xlabel "Cycle" font ",16"
set key font ",16" outside right

# Title for data index i: pulls the matching "# p = …, μ = …" header line
# from the file (data index i corresponds to the (i+1)-th such header line),
# keeping only the "p = …" field.
blocktitle(datafile, i) = system(sprintf("grep -a '^# p =' '%s' | sed -n '%dp' | sed 's/^# //; s/,.*//'", datafile, i + 1))

nblocks(datafile) = system(sprintf("grep -ac '^# p =' '%s'", datafile)) + 0

# == Number of Pauli strings ============================================================
set output "output/Resources_prob_tm_NP_pp_graph.png"
datafile = "output/carlos_varseed_pp/resources_prob_tm_NP_pp.txt"
set title "46-qubit XXZ circuit — PP Pauli string number" font ",20"
set ylabel "N. of Pauli strings" font ",16"
set logscale y
nb = nblocks(datafile)
plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title blocktitle(datafile, i)

# == Mean Pauli weight ==================================================================
set output "output/Resources_prob_tm_PW_pp_graph.png"
datafile = "output/carlos_varseed_pp/resources_prob_tm_PW_pp.txt"
set title "46-qubit XXZ circuit — PP mean Pauli weight" font ",20"
set ylabel "Mean avg. Pauli weight" font ",16"
set logscale xy
set yrange [.9:10]
set xrange [1:25]
nb = nblocks(datafile)
plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title blocktitle(datafile, i)

# == Mean |Pauli coefficient| ===========================================================
set output "output/Resources_prob_tm_PC_pp_graph.png"
datafile = "output/carlos_varseed_pp/resources_prob_tm_PC_pp.txt"
set title "46-qubit XXZ circuit — PP mean Pauli coeff." font ",20"
set ylabel "Mean avg. |Pauli coeff.|" font ",16"
unset xrange
unset yrange
unset logscale xy
set logscale y
nb = nblocks(datafile)
plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title blocktitle(datafile, i)
