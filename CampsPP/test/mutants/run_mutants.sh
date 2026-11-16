#!/bin/bash
# T10 mutant runner: one fresh `Pkg.test()` per mutant, on a scratch copy of the repo.
#
# usage: bash run_mutants.sh <scratch dir> [cap seconds] [ID ...]
#   <scratch dir>  where the copy (work/) and the logs (results/) go; anything outside the repo,
#                  e.g. "$TMPDIR/mutants". work/ is recreated from the repo on every call.
#   cap seconds    wall-clock cap per run (default 600); a run over the cap is killed together
#                  with Pkg.test's child julia and recorded as TIMEOUT.
#   ID ...         mutants to run (patches/<ID>.patch); "BASE" runs the unmutated copy.
#                  Default: BASE, then every patch in patches/.
#
# For each mutant: apply the patch to the copy, run Pkg.test() in the copy's environment (a fresh
# process, so the mutated source is precompiled anew), revert, and check that the copy's src equals
# the repo's again (the runner stops otherwise). The repo itself is never modified.
# Run one at a time: each run is one Julia test process (plus Pkg.test's parent).
#
# Results: results/summary.tsv (id, status, pass, fail, error, broken, seconds) and one log per
# run. Status: killed (failures or errors), SURVIVED (suite passed), TIMEOUT (cap hit; the log
# still has the failures reported before the hang), CRASH (no test summary), clean (BASE passed).
# who_caught.py <results dir> ID ... lists which tests caught each mutant.
set -m # background jobs get their own process group, so the cap can also kill Pkg.test's child
HERE=$(cd "$(dirname "$0")" && pwd -P)
REPO=$(cd "$HERE/../../.." && pwd -P)
PD=$HERE/patches
[[ $# -ge 1 ]] || { sed -n '2,20p' "$0"; exit 2; }
mkdir -p "$1"; S=$(cd "$1" && pwd -P); shift
CAP=${1:-600}; [[ $# -ge 1 ]] && shift
[[ $CAP =~ ^[0-9]+$ ]] || { echo "cap must be a number of seconds, got '$CAP'"; exit 2; }
if [[ $# -ge 1 ]]; then IDS=("$@"); else
  IDS=(BASE); for p in "$PD"/*.patch; do IDS+=("$(basename "$p" .patch)"); done
fi
case $S/ in "$REPO"/*) echo "scratch dir must be outside the repo"; exit 2;; esac

# Pkg.test builds its test environment under TMPDIR and recomputes the relative [sources] path
# to dmrgCliffordMPS from there; a symlinked TMPDIR (macOS /tmp -> /private/tmp) or a symlinked
# dmrgCliffordMPS makes that check fail (AssertionError on normpath). Use physical paths and copies.
export TMPDIR=$(cd "${TMPDIR:-/tmp}" && pwd -P)
W=$S/work; OUT=$S/results; SUM=$OUT/summary.tsv
rm -rf "$W"; mkdir -p "$W" "$OUT"
cp -R "$REPO/CampsPP" "$W/CampsPP"
cp -R "$REPO/dmrgCliffordMPS" "$W/dmrgCliffordMPS"
[[ -f $SUM ]] || printf 'id\tstatus\tpass\tfail\terror\tbroken\tseconds\n' > "$SUM"

clean_or_die() {
  if ! diff -r -q "$W/CampsPP/src" "$REPO/CampsPP/src" >/dev/null; then
    echo "runner: copy src differs from the repo after $1; stopping"; exit 1
  fi
}

clean_or_die "setup"
for id in "${IDS[@]}"; do
  log=$OUT/$id.log
  if [[ $id != BASE ]]; then
    patch -p1 -s -d "$W" -i "$PD/$id.patch" || { echo "runner: patch $id does not apply; regenerate with gen_patches.py"; exit 1; }
  fi
  echo "start $id $(date +%T)"
  t0=$(date +%s)
  julia --project="$W/CampsPP" -e 'using Pkg; Pkg.test()' > "$log" 2>&1 &
  pid=$!
  status=""
  while kill -0 $pid 2>/dev/null; do
    sleep 5
    if (( $(date +%s) - t0 > CAP )); then
      kill -TERM -- -$pid 2>/dev/null; sleep 5; kill -KILL -- -$pid 2>/dev/null
      status=TIMEOUT; break
    fi
  done
  wait $pid 2>/dev/null; rc=$?
  secs=$(( $(date +%s) - t0 ))
  [[ $id != BASE ]] && patch -p1 -R -s -d "$W" -i "$PD/$id.patch"
  clean_or_die "$id"
  counts=$(grep -oE 'Some tests did not pass: [0-9]+ passed, [0-9]+ failed, [0-9]+ errored, [0-9]+ broken' "$log" | head -1 |
           sed -E 's/.*: ([0-9]+) passed, ([0-9]+) failed, ([0-9]+) errored, ([0-9]+) broken/\1\t\2\t\3\t\4/')
  if [[ -z $status ]]; then
    if [[ -n $counts ]]; then status=killed
    elif grep -q 'tests passed' "$log"; then status=$([[ $id == BASE ]] && echo clean || echo SURVIVED)
    else status=CRASH; fi
  fi
  if [[ -z $counts ]]; then
    p=$(grep -E '^CampsPP package +\|' "$log" | tail -1 | awk -F'|' '{print $2}' | awk '{print $1}')
    nf=$(grep -cE ': (Test Failed|Error During Test) at' "$log")
    counts="${p:--1}\t${nf}\t-\t-"
  fi
  printf "%s\t%s\t$counts\t%s\n" "$id" "$status" "$secs" >> "$SUM"
  echo "$id: $status ($secs s, rc=$rc) $(echo -e "$counts" | tr '\t' ' ')"
done
echo "runner: done"
