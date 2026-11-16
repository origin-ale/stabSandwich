#!/bin/bash
# T11 sensitivity check: replay a few Phase 0 CAMPS-PP samples with a known defect and
# confirm that compare_prob_tm.jl flags it through the ⟨O⟩/switch checks (not only
# through the input check, which a moved doping mask trips trivially). A clean replay of
# the same samples must pass. Run locally from the repository root:
#
#   scripts/t11/sensitivity.sh <Phase 0 CAMPS-PP run dir> <scratch dir>
#
# Defects (replay_prob_tm.jl `inject=`): `phase`, the magic phase × 1.05; `mask`, every
# doped XX-YY pair moved one XXZ triple later. p = 0 can't detect either (nothing is
# doped), so the samples are switched ones at p = 0.02 and 0.04 plus an unswitched one
# at p = 0.01 (CAMPS all the way, 1e-10 tolerance).
set -uo pipefail

run=${1:?usage: sensitivity.sh <Phase 0 run dir> <scratch dir>}
scratch=${2:?usage: sensitivity.sh <Phase 0 run dir> <scratch dir>}
p=${T11_SENS_P:-0.01,0.02,0.04}
samples=${T11_SENS_SAMPLES:-1,6,7}
old="$run/files/output"
fail=0

for inject in none phase mask; do
  out="$scratch/sens_$inject"
  echo "== inject=$inject (p=$p, samples=$samples)"
  julia --project=. scripts/t11/replay_prob_tm.jl "run=$run" methods=campspp "p=$p" \
    "samples=$samples" "inject=$inject" "out=$out" >"$out.replay.log" 2>&1 \
    || { echo "replay failed, see $out.replay.log"; fail=1; continue; }
  julia --project=. scripts/t11/compare_prob_tm.jl "old=$old" "new=$out" >"$out.compare.log" 2>&1
  code=$?
  # Violations other than the input check
  nout=$(grep '^VIOLATION' "$out.compare.log" | grep -cv -e 'doped gates differ' -e 'initial states differ')
  nin=$(grep -c -e '^VIOLATION.*doped gates differ' "$out.compare.log")
  echo "   compare exit $code; output violations $nout; input violations $nin"
  if [ "$inject" = none ]; then
    [ "$code" -eq 0 ] || { echo "   FAIL: the clean replay must pass"; fail=1; }
  else
    [ "$code" -eq 1 ] && [ "$nout" -gt 0 ] || { echo "   FAIL: defect not flagged by the output checks"; fail=1; }
    # Per (p, sample): which cases the output checks flagged
    grep '^VIOLATION' "$out.compare.log" | grep -v 'differ$' \
      | sed -E 's/^VIOLATION: (campspp p = [0-9.]+, sample [0-9]+).*/\1/' | sort | uniq -c | sed 's/^/   /'
  fi
done
[ "$fail" -eq 0 ] && echo "Sensitivity check PASSED" || echo "Sensitivity check FAILED"
exit "$fail"
