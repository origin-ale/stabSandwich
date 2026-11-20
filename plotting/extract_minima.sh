#!/bin/sh
# Extract, from each M*_switch_optimization_avgs_<N>.txt file, the row with the
# smallest computation time, and collect them (prefixed by M) in switch_minima.txt
#
#   plotting/extract_minima.sh [DIR [N]]
#
# DIR defaults to output/switch_optimization (where switch_optimization.jl writes;
# old data: output/old/carlos_switch), N to 20.

dir="${1:-output/switch_optimization}"
N="${2:-20}"
out="$dir/switch_minima.txt"

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
