#!/bin/bash
# Switch-optimization times with the compile time taken out, for the plot scripts
# (plot_switches.plt, plot_switches_box.plt). Run from the repo root.
#
#   switch_times.sh samples FULL   rows "c t_1 … t_n", as in FULL
#   switch_times.sh avgs FULL      rows "c mean std.err", as in the _avgs_ file
#                                  (std.err = std/sqrt(n), std with n-1 like Julia's)
#   switch_times.sh hascompile FULL   prints 1 if FULL and its compile file both exist, else 0
#
# FULL is an M<M>_switch_optimization_full_<N>.txt file. If the matching
# M<M>_switch_optimization_compile_<N>.txt exists (switch_optimization.jl writes it
# in the same layout), each sample has its compile time subtracted; otherwise the
# samples are used as they are. A switch point the compile file has no row for (a
# run stopped between writing the two files) is left out, with a warning.
set -u
die() { echo "switch_times.sh: $*" >&2; exit 1; }
[ $# -eq 2 ] || die "usage: switch_times.sh samples|avgs|hascompile FULL"
mode=$1 full=$2
case "$full" in *_switch_optimization_full_*) ;; *) die "not a _full_ file: $full" ;; esac
dir=$(dirname "$full")
compile="$dir/$(basename "$full" | sed 's/_switch_optimization_full_/_switch_optimization_compile_/')"

case "$mode" in
hascompile) [ -f "$full" ] && [ -f "$compile" ] && echo 1 || echo 0; exit 0 ;;
samples|avgs) ;;
*) die "unknown mode $mode" ;;
esac
[ -f "$full" ] || die "no $full"

samples() {
  if [ -f "$compile" ]; then
    awk -v cf="$compile" -v ff="$full" '
      !/^#/ && NF > 1 {
        if (FILENAME == cf) { nc[$1] = NF; for (i = 2; i <= NF; i++) comp[$1, i] = $i; next }
        if (!($1 in nc)) { printf "switch_times.sh: no compile row for c = %s in %s, left out\n", $1, cf > "/dev/stderr"; next }
        if (nc[$1] != NF) { printf "switch_times.sh: c = %s has %d samples in %s but %d in %s\n", $1, NF - 1, ff, nc[$1] - 1, cf > "/dev/stderr"; bad = 1; exit 1 }
        line = $1
        for (i = 2; i <= NF; i++) line = line " " sprintf("%.9g", $i - comp[$1, i])
        print line
      }
      END { if (bad) exit 1 }' "$compile" "$full"
  else
    awk '!/^#/ && NF > 1' "$full"
  fi
}

if [ "$mode" = samples ]; then
  samples
else
  samples | awk '{ n = NF - 1; s = 0; for (i = 2; i <= NF; i++) s += $i; m = s / n
                   v = 0; for (i = 2; i <= NF; i++) v += ($i - m)^2
                   printf "%s %.9g %.9g\n", $1, m, (n > 1 ? sqrt(v / (n - 1) / n) : 0) }'
fi
