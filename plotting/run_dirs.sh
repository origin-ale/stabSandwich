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
#   run_dirs.sh method FILE METHOD
#       Prints METHOD's 0-based position in the "methods=…" list of FILE's
#       first (header) line (methods=a,b or Julia's methods=[:a, :b]), -1 if
#       the list doesn't include it. A FILE without the list gives 0.
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
method)
  [ $# -eq 3 ] || die "usage: run_dirs.sh method FILE METHOD"
  [ -f "$2" ] || die "no $2"
  # methods=campspp,pp, or Julia's methods=[:campspp, :pp] in older files
  ms=$(head -1 "$2" | grep -oE 'methods=(\[[^]]*\]|[^ ]*)' | cut -d= -f2 | tr -d '[]: ')
  [ -n "$ms" ] || { echo 0; exit 0; }
  echo "$ms" | tr ',' '\n' | awk -v m="$3" '$0 == m { p = NR - 1 } END { print (p == "" ? -1 : p) }'
  ;;
*) die "usage: run_dirs.sh data DATADIR [RUN] | plotdir DIR... | method FILE METHOD" ;;
esac
