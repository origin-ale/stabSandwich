#!/bin/sh
# Extract, from each M*_switch_optimization_avgs_20.txt file, the row with the
# smallest computation time, and collect them (prefixed by M) in switch_minima.txt

dir="${1:-output/switch}"
out="$dir/switch_minima.txt"

{
    echo "# Minima of computation time"
    echo "# M c Time (s) Error"
    for f in "$dir"/M*_switch_optimization_avgs_20.txt; do
        M=$(basename "$f" | sed 's/^M\([0-9]*\)_.*/\1/')
        awk -v M="$M" '
            /^[ \t]*#/ || NF < 3 { next }
            best == "" || $2 < best { best = $2; row = $0 }
            END { if (row != "") print M, row }
        ' "$f"
    done | sort -n -k1
} > "$out"

echo "Wrote $out"
