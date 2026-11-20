#!/bin/bash
# Run every Comptime_prob_tm plot script on one fetched Comptime_prob_tm run.
# Run from the repo root:
#
#   plotting/plot_prob_tm_all.sh            # the newest fetched Comptime_prob_tm run
#   plotting/plot_prob_tm_all.sh RUN-ID     # a given one
#
# Plots go to output/plots/RUN-ID/. The run's own settings are passed on:
# each method's position in its `methods` (j_campspp, j_pp) and the threshold
# suffix of the evs_prob_tm_<method>_<suffix>.txt files. Scripts for a method
# the run didn't include are skipped. A script that fails doesn't stop the
# others; the exit status is non-zero if any failed.
set -uo pipefail

die() { echo "plot_prob_tm_all.sh: $*" >&2; exit 1; }
[ -d plotting ] && [ -d output/runs ] || die "run from the repo root (no plotting/ or output/runs/ here)"
sub=output/Comptime_prob_tm

if [ $# -ge 1 ]; then
  id=$1
else
  # Runs before the per-script output folders wrote to output/ directly, so
  # require the folder, not just the name.
  id=""
  for r in $(ls -1 output/runs | grep -- '-Comptime_prob_tm$' | LC_ALL=C sort -r); do
    [ -d "output/runs/$r/files/$sub" ] && { id=$r; break; }
  done
  [ -n "$id" ] || die "no fetched Comptime_prob_tm run with $sub in output/runs"
fi
dir="output/runs/$id/files/$sub"
[ -d "$dir" ] || die "run $id has no $sub (no $dir)"
echo "Run: $id"

# Method positions, from the "# methods=… " header of comptime_prob_tm.txt
methods=$(head -1 "$dir/comptime_prob_tm.txt" 2>/dev/null | grep -o 'methods=[^ ]*' | cut -d= -f2)
[ -n "$methods" ] || die "no methods= header in $dir/comptime_prob_tm.txt"
pos() { echo "$methods" | tr ',' '\n' | grep -nx "$1" | cut -d: -f1; }
jc=$(pos campspp); jp=$(pos pp)
echo "Methods: $methods"

# Threshold suffix of the evs files (one per method expected)
thl() {
  local s
  s=$(ls "$dir" | sed -n "s/^evs_prob_tm_$1_\(.*\)\.txt\$/\1/p" | grep -v '_full$')
  [ "$(echo "$s" | grep -c .)" -eq 1 ] && echo "$s"
}

ran=(); skipped=(); failed=()
plot() {   # plot NAME NEEDS [gnuplot -e assignments]
  local name=$1 needs=$2 extra=${3:-}
  local m
  for m in $needs; do
    case ",$methods," in *",$m,"*) ;; *) skipped+=("$name (no $m)"); return ;; esac
  done
  echo "--- $name"
  if gnuplot -e "run='$id'; $extra" "plotting/$name.plt"; then ran+=("$name"); else failed+=("$name"); fi
}

plot boxplot_resources_prob_tm_campspp campspp
plot plot_resources_prob_tm_campspp    campspp
plot boxplot_resources_prob_tm_pp      pp
plot plot_resources_prob_tm_pp         pp
plot plot_resources_peaks_pp           pp
plot plot_comptime_peaks_pp            pp
if [ -n "$jc" ] && [ -n "$jp" ]; then
  both="j_campspp=$((jc - 1)); j_pp=$((jp - 1))"
  plot plot_comptime_prob_tm     "campspp pp" "$both"
  plot plot_comptime_prob_tm_box "campspp pp" "$both"
  tc=$(thl campspp); tp=$(thl pp)
  if [ -n "$tc" ] && [ -n "$tp" ]; then
    plot plot_evs_diff_p050 "campspp pp" "thl_campspp='$tc'; thl_pp='$tp'"
  else
    skipped+=("plot_evs_diff_p050 (not exactly one evs_prob_tm_{campspp,pp}_<thl>.txt each)")
  fi
else
  skipped+=("plot_comptime_prob_tm, plot_comptime_prob_tm_box, plot_evs_diff_p050 (need both campspp and pp)")
fi

echo
echo "Plots in output/plots/$id/"
echo "  ran:     ${#ran[@]}"
[ ${#skipped[@]} -eq 0 ] || printf '  skipped: %s\n' "${skipped[@]}"
[ ${#failed[@]} -eq 0 ] || { printf '  FAILED:  %s\n' "${failed[@]}"; exit 1; }
