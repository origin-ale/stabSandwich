set term png size 850, 600 font ",16"

set xlabel "Cycle" font ",16"
set key font ",16" outside right

# Title for data index i: pulls the matching "# p = …, μ = …" header line
# from the file (data index i corresponds to the (i+1)-th such header line),
# keeping only the "p = …" field.
blocktitle(datafile, i) = system(sprintf("grep -a '^# p =' '%s' | sed -n '%dp' | sed 's/^# //; s/,.*//'", datafile, i + 1))

nblocks(datafile) = system(sprintf("grep -ac '^# p =' '%s'", datafile)) + 0

# Count of in-range rows in block `blk` (0-indexed) with a non-zero error bar;
# single-sample cycles have err = 0 and would blow up a weighted (yerrors) fit,
# so they're excluded from weighting and must be counted separately, since
# gnuplot's stats aborts outright (rather than reporting zero) when every row
# in a NaN-filtered column is invalid.
nzcount(datafile, blk, lo, hi) = system(sprintf("awk 'BEGIN{b=-1} /^# p =/{b++} b==%d && $1>=%d && $1<=%d && $3>0 {c++} END{print c+0}' '%s'", blk, lo, hi, datafile)) + 0

set fit quiet nolog errorvariables

# == Number of Pauli strings ============================================================
set output "output/Resources_prob_tm_NP_campspp_graph.png"
datafile = "output/carlos_varseed64/resources_prob_tm_NP_campspp.txt"
set title "46-qubit XXZ circuit — CAMPS-PP Pauli string number" font ",20"
set ylabel "N. of Pauli strings" font ",16"
set logscale y
nb = nblocks(datafile)
plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title blocktitle(datafile, i)

# == Mean Pauli weight ==================================================================
set output "output/Resources_prob_tm_PW_campspp_graph.png"
datafile = "output/carlos_varseed64/resources_prob_tm_PW_campspp.txt"
set title "46-qubit XXZ circuit — CAMPS-PP mean Pauli weight" font ",20"
set ylabel "Mean avg. Pauli weight" font ",16"
set logscale xy
set yrange [.9:10]
set xrange [9:25]
nb = nblocks(datafile)

# Power-law fit w(c) = a·c^γ per p block, over the plotted cycle range only.
fitmin = 17
fitmax = 25
f(x) = a*x**gamma
array pw_a[nb]
array pw_gamma[nb]
print sprintf("\n-- %s --", datafile)
do for [i=1:nb] {
    # Low-p blocks only survive a few cycles before the window opens, so they
    # have < 3 points in range (or are still flat); γ is undetermined there.
    # Mark them NaN so both the curve and its key entry drop out below. Also
    # require >= 3 non-zero-error points, since those are all the weighted
    # fit can use.
    stats [fitmin:fitmax] datafile index i-1 using 1:2 nooutput
    nz = nzcount(datafile, i-1, fitmin, fitmax)
    if (STATS_records >= 3 && STATS_max_y > STATS_min_y && nz >= 3) {
        a = 1.0
        gamma = 0.5
        fit [fitmin:fitmax] f(x) datafile index i-1 using 1:2:($3 > 0 ? $3 : NaN) yerrors via a, gamma
        pw_a[i] = a
        pw_gamma[i] = gamma
        print sprintf("  %-14s a = %.4f +/- %.4f   gamma = %.4f +/- %.4f   chi2/ndf = %.3f (ndf=%d)", \
            blocktitle(datafile, i-1), a, a_err, gamma, gamma_err, FIT_WSSR/FIT_NDF, FIT_NDF)
    } else {
        pw_a[i] = NaN
        pw_gamma[i] = NaN
        print sprintf("  %-14s skipped (insufficient in-range data)", blocktitle(datafile, i-1))
    }
}

# γ rides along in the data series' own key entry, so the exponent sits next to
# the p it belongs to and the fit curves add no entries of their own.
fitlabel(i) = (pw_gamma[i+1] == pw_gamma[i+1]) \
    ? sprintf("%s, γ = %.2f", blocktitle(datafile, i), pw_gamma[i+1]) \
    : blocktitle(datafile, i)

plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title fitlabel(i), \
     for [i=0:nb-1] (x >= fitmin && x <= fitmax) ? pw_a[i+1]*x**pw_gamma[i+1] : NaN \
         with lines dt 2 lw 2 lc rgb "gray" notitle

# == Mean |Pauli coefficient| ===========================================================
set output "output/Resources_prob_tm_PC_campspp_graph.png"
datafile = "output/carlos_varseed64/resources_prob_tm_PC_campspp.txt"
set title "46-qubit XXZ circuit — CAMPS-PP mean Pauli coeff." font ",20"
set ylabel "Mean avg. |Pauli coeff.|" font ",16"
unset xrange
unset yrange
unset logscale xy
set logscale y
nb = nblocks(datafile)
plot for [i=0:nb-1] datafile index i using 1:2:3 with yerrorlines lc i+1 title blocktitle(datafile, i)
