#!/bin/sh
# Extract, from each M*_switch_optimization_avgs_<N>.txt file, the row with the
# smallest computation time, and collect them (prefixed by M) in switch_minima.txt
#
#   plotting/extract_minima.sh [DIR [N]]
#
# DIR defaults to the newest fetched run's copy of output/switch_optimization
# (output/switch_optimization itself if no run was fetched; see run_dirs.sh;
# old data: output/old/carlos_switch), N to 20.
#
# switch_minima.txt goes to DIR's plot folder, output/plots/<run id>, where
# plot_minima.plt and plot_minima_c.plt read it: derived data stays out of the
# run's checksummed copy.

# Run from the repo root, like the plot scripts (run_dirs.sh paths are relative to it).
rd="$(dirname "$0")/run_dirs.sh"
dir="${1:-$("$rd" data output/switch_optimization)}" || exit 1
N="${2:-20}"
plotdir=$("$rd" plotdir "$dir") || exit 1
out="$plotdir/switch_minima.txt"

{
    echo "# Minima of computation time"
    echo "# M c Time (s) Error"
    for f in "$dir"/M*_switch_optimization_avgs_"$N".txt; do
        [ -f "$f" ] || continue
        M=$(basename "$f" | sed 's/^M\([0-9]*\)_.*/\1/')
        awk -v M="$M" '
            /^[ \t]*#/ || NF < 3 { next }
            best == "" || $2 < best { best = $2; row = $0 }
            END { if (row != "") print M, row }
        ' "$f"
    done | sort -n -k1
} > "$out"

echo "Wrote $out"
