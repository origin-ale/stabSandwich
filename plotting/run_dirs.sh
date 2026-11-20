#!/bin/bash
# Run selection for the plot scripts: which data to read, and where the plots
# go. Called by plotting/run.plt (the gnuplot scripts), extract_minima.sh and
# plot_tm_vp_v2.jl. Run from the repo root, like the plot scripts.
#
#   run_dirs.sh data DATADIR [RUN]
#       Input directory for data a simulation script writes to DATADIR (e.g.
#       output/Comptime_prob_tm): the run's own copy, output/runs/RUN/files/DATADIR,
#       which fetch.sh installs and checks against the run's SHA256SUMS. Without
#       RUN, the newest fetched run that wrote DATADIR; if no fetched run did,
#       DATADIR itself (data not fetched through fetch.sh).
#   run_dirs.sh plotdir DIR...
#       Creates and prints output/plots/<id>, the plot folder for the input
#       DIR(s): <id> is the run id for a DIR under output/runs/<id>/, otherwise
#       the DIR's own name (output/old/carlos_varseed -> carlos_varseed). DIRs
#       from different runs give "<id1>+<id2>".
set -u
die() { echo "run_dirs.sh: $*" >&2; exit 1; }
norm() { local d=${1#./}; echo "${d%/}"; }

case "${1:-}" in
data)
  { [ $# -ge 2 ] && [ $# -le 3 ]; } || die "usage: run_dirs.sh data DATADIR [RUN]"
  datadir=$(norm "$2")
  if [ $# -eq 3 ]; then
    d="output/runs/$3/files/$datadir"
    [ -d "$d" ] || die "run $3 has no $datadir (no $d)"
    echo "$d"; exit 0
  fi
  for r in $(ls -1 output/runs 2>/dev/null | LC_ALL=C sort -r); do
    d="output/runs/$r/files/$datadir"
    [ -d "$d" ] && { echo "$d"; exit 0; }
  done
  echo "$datadir"
  ;;
plotdir)
  shift
  [ $# -ge 1 ] || die "usage: run_dirs.sh plotdir DIR..."
  id=""
  for d in "$@"; do
    d=$(norm "$d")
    case "$d" in
      output/runs/*/*) i=${d#output/runs/}; i=${i%%/*} ;;
      *) i=$(basename "$d") ;;
    esac
    case "+$id+" in *"+$i+"*) ;; *) id=${id:+$id+}$i ;; esac
  done
  out="output/plots/$id"
  mkdir -p "$out" || die "cannot create $out"
  echo "plots -> $out" >&2
  echo "$out"
  ;;
*) die "usage: run_dirs.sh data DATADIR [RUN] | plotdir DIR..." ;;
esac
