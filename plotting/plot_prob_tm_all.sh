#!/bin/bash
# Run every Comptime_prob_tm plot script on one fetched Comptime_prob_tm run.
# Run from the repo root:
#
#   plotting/plot_prob_tm_all.sh            # the newest fetched Comptime_prob_tm run
#   plotting/plot_prob_tm_all.sh RUN-ID     # a given one
#
# Plots go to output/plots/RUN-ID/. The threshold suffix of the run's
# evs_prob_tm_<method>_<suffix>.txt files is passed on; the computation-time
# scripts find each method's column from the run's `methods` themselves and
# plot whichever of CAMPS-PP and PP the run has. Scripts for a method the run
# didn't include are skipped. A script that fails doesn't stop the others; the
# exit status is non-zero if any failed.
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

# The run's methods, from the "# methods=… " header of comptime_prob_tm.txt
methods=$(head -1 "$dir/comptime_prob_tm.txt" 2>/dev/null | grep -o 'methods=[^ ]*' | cut -d= -f2)
[ -n "$methods" ] || die "no methods= header in $dir/comptime_prob_tm.txt"
has() { case ",$methods," in *",$1,"*) return 0 ;; esac; return 1; }
echo "Methods: $methods"

# Threshold suffix of the evs files (one per method expected)
thl() {
  local s
  s=$(ls "$dir" | sed -n "s/^evs_prob_tm_$1_\(.*\)\.txt\$/\1/p" | grep -v '_full$')
  [ "$(echo "$s" | grep -c .)" -eq 1 ] && echo "$s"
}

ran=(); skipped=(); failed=()
plot() {   # plot NAME [gnuplot -e assignments]
  echo "--- $1"
  if gnuplot -e "run='$id'; ${2:-}" "plotting/$1.plt"; then ran+=("$1"); else failed+=("$1"); fi
}

if has campspp; then
  plot boxplot_resources_prob_tm_campspp
  plot plot_resources_prob_tm_campspp
else
  skipped+=("boxplot_resources_prob_tm_campspp, plot_resources_prob_tm_campspp (no campspp)")
fi
if has pp; then
  plot boxplot_resources_prob_tm_pp
  plot plot_resources_prob_tm_pp
  plot plot_resources_peaks_pp
  plot plot_comptime_peaks_pp
else
  skipped+=("boxplot_resources_prob_tm_pp, plot_resources_prob_tm_pp, plot_resources_peaks_pp, plot_comptime_peaks_pp (no pp)")
fi
# Computation time: whichever of the two methods the run has
if has campspp || has pp; then
  plot plot_comptime_prob_tm
  plot plot_comptime_prob_tm_box
else
  skipped+=("plot_comptime_prob_tm, plot_comptime_prob_tm_box (neither campspp nor pp)")
fi
# Expectation values: a comparison, so both
if has campspp && has pp; then
  tc=$(thl campspp); tp=$(thl pp)
  if [ -n "$tc" ] && [ -n "$tp" ]; then
    plot plot_evs_diff_p050 "thl_campspp='$tc'; thl_pp='$tp'"
  else
    skipped+=("plot_evs_diff_p050 (not exactly one evs_prob_tm_{campspp,pp}_<thl>.txt each)")
  fi
else
  skipped+=("plot_evs_diff_p050 (compares campspp with pp, needs both)")
fi

echo
echo "Plots in output/plots/$id/"
echo "  ran:     ${#ran[@]}"
[ ${#skipped[@]} -eq 0 ] || printf '  skipped: %s\n' "${skipped[@]}"
[ ${#failed[@]} -eq 0 ] || { printf '  FAILED:  %s\n' "${failed[@]}"; exit 1; }
