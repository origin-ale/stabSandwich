#!/bin/bash
# Launch a Julia script on the remote machine, detached, one run at a time.
# Run this ON THE TARGET, from anywhere inside the repo clone.
#
#   scripts/remote/run.sh [-t THREADS] [--no-pull] SCRIPT [ARGS...]
#   scripts/remote/run.sh --status     # current run (if any) and its log tail
#   scripts/remote/run.sh --list       # all runs with their status
#
# What a run does:
#   1. takes the lock output/runs/.lock (refuses if another run is alive),
#   2. refuses if tracked files have local changes (the commit must describe the code),
#   3. git pull --ff-only (skip with --no-pull),
#   4. starts `julia --project=. -t THREADS SCRIPT ARGS...` under nohup,
#   5. when Julia exits: lists every file under output/ that the run created or
#      modified (files.txt), copies them into the run directory, finalizes the
#      run record, writes SHA256SUMS, and releases the lock.
#
# Everything about a run lives in output/runs/<run-id>/:
#   record.txt  commit, host, threads, command, start/end time, exit code, status
#   log.txt     Julia's stdout and stderr
#   files.txt   output files written by the run (paths relative to the repo root)
#   files/      a copy of those files at the same relative paths. Scripts write to
#               fixed names under output/, so the next run overwrites them; this
#               copy keeps every run fetchable until you delete it.
#   SHA256SUMS  checksums (paths relative to the run directory) that fetch.sh
#               uses to reject truncated transfers; written last, so its
#               presence means the run has finished
#
# THREADS defaults to $JULIA_NUM_THREADS, else 1. It must be a number, or N,M
# (default and interactive pools). "auto" is refused, because it means every core.
# Every script gives BLAS and Strided 1 thread unless a `blas_threads=K` argument
# says otherwise (plan item 10, 2026-11-20); a K above 1 adds K BLAS threads, so
# THREADS plus K must stay within 20, the usage limit on the target.
# Comptime_prob_tm.jl and the T11 replay run samples one after another and need no
# more than -t 1 (more Julia threads only slow them down); give
# TransferredMag_varparams.jl, which runs one sample per Julia thread, a larger -t.
# The julia binary is $JULIA if set, else `julia` on PATH.
#
# Stop a run with `kill <worker pid>` (shown by --status); Julia gets SIGTERM,
# then SIGKILL if it is still running 10 s later. The run is recorded as killed
# and the lock is released.
#
# Runs are attributed by modification time, which is sound only because the
# lock keeps runs from overlapping. Don't write into output/ by hand during a run.
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
self="$repo/scripts/remote/run.sh"
runs="$repo/output/runs"
lock="$runs/.lock"

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }
die() { echo "run.sh: $*" >&2; exit 1; }

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"; else shasum -a 256 "$@"; fi
}

ncores() { nproc 2>/dev/null || getconf _NPROCESSORS_ONLN 2>/dev/null || echo "?"; }

field() { sed -n "s/^$2: //p" "$1" | tail -1; }   # field FILE KEY

max_threads=20       # usage limit on the target (user, 2026-11-16)
default_threads=1    # the serial scripts are fastest at -t 1 (plan item 10, 2026-11-20)

total_threads() {    # "8" -> 8, "8,1" -> 9; fails on anything else ("auto", "")
  [[ $1 =~ ^[0-9]+(,[0-9]+)?$ ]] || return 1
  local a=${1%%,*} b=0
  [[ $1 == *,* ]] && b=${1#*,}
  echo $((10#$a + 10#$b))
}

lock_owner() {   # prints "pid run-id" of a live lock holder, or nothing
  [ -d "$lock" ] || return 0
  local pid id
  pid=$(cat "$lock/pid" 2>/dev/null) || return 0
  id=$(cat "$lock/run" 2>/dev/null)
  if [ -n "$pid" ] && kill -0 "$pid" 2>/dev/null; then echo "$pid $id"; fi
}

# ---------------------------------------------------------------- worker ---
# Re-invoked by the launcher under nohup; owns the lock until it exits.
if [ "${1:-}" = "__worker" ]; then
  rundir=$2; threads=$3; shift 3
  cd "$repo" || exit 1
  record="$rundir/record.txt"
  jpid=""; status=""; reaper=""
  grace=10   # seconds Julia gets to exit on SIGTERM before it gets SIGKILL

  # Julia can hang after SIGTERM (it prints its backtrace, then never exits, e.g.
  # with several threads busy), so a timer sends SIGKILL after $grace seconds.
  # The subshell does not inherit the EXIT trap, so it never touches the lock.
  on_signal() {
    status="killed"
    [ -n "$jpid" ] || return 0
    kill -TERM "$jpid" 2>/dev/null
    [ -n "$reaper" ] && return 0
    ( sleep "$grace"
      kill -0 "$jpid" 2>/dev/null || exit 0
      echo "julia still running $grace s after SIGTERM; sent SIGKILL" >>"$record"
      kill -KILL "$jpid" 2>/dev/null ) &
    reaper=$!
  }
  trap on_signal TERM INT
  trap 'rm -rf "$lock"' EXIT
  echo $$ >"$lock/pid"
  echo "worker pid: $$" >>"$record"

  touch "$rundir/.start"
  sleep 1   # keep files written in the same second as the marker unambiguous

  "${JULIA:-julia}" --project=. -t "$threads" "$@" >"$rundir/log.txt" 2>&1 </dev/null &
  jpid=$!
  echo "julia pid: $jpid" >>"$record"
  wait "$jpid"; code=$?
  # A trapped signal interrupts `wait`; wait again for Julia to actually exit.
  while kill -0 "$jpid" 2>/dev/null; do wait "$jpid"; code=$?; done
  # Julia is reaped: stop the timer before its pid check could hit a reused pid.
  [ -z "$reaper" ] || kill "$reaper" 2>/dev/null
  [ -n "$status" ] || { [ "$code" -eq 0 ] && status="done" || status="failed"; }

  (cd "$repo" && find output -path output/runs -prune -o -type f ! -name .DS_Store \
     -newer "$rundir/.start" -print | LC_ALL=C sort) >"$rundir/files.txt"
  others=$(cd "$repo" && find . \( -path ./.git -o -path ./output \) -prune -o -type f \
     -newer "$rundir/.start" -print | LC_ALL=C sort | tr '\n' ' ')
  rm -f "$rundir/.start"

  while IFS= read -r f; do
    mkdir -p "$rundir/files/$(dirname "$f")" && cp -p "$f" "$rundir/files/$f" \
      || echo "warning: could not copy $f into the run directory" >>"$record"
  done <"$rundir/files.txt"

  {
    echo "end: $(date '+%Y-%m-%d %H:%M:%S %z')"
    echo "exit code: $code"
    echo "output files: $(wc -l <"$rundir/files.txt" | tr -d ' ')"
    [ -z "$others" ] || echo "warning: files modified outside output/ (not fetched): $others"
    echo "status: $status"
  } >>"$record"

  cd "$rundir" || exit 1
  { sed 's|^|files/|' files.txt; printf '%s\n' log.txt record.txt files.txt; } \
    | while IFS= read -r f; do sha256 "$f"; done >SHA256SUMS.tmp
  mv SHA256SUMS.tmp SHA256SUMS
  exit 0
fi

# ---------------------------------------------------------- status / list ---
if [ "${1:-}" = "--list" ]; then
  [ -d "$runs" ] || exit 0
  owner=$(lock_owner); live=${owner#* }
  for d in $(ls -1 "$runs" 2>/dev/null | LC_ALL=C sort); do
    r="$runs/$d/record.txt"; [ -f "$r" ] || continue
    s=$(field "$r" status)
    if [ -z "$s" ]; then [ "$d" = "$live" ] && s="running" || s="interrupted"; fi
    printf '%-40s %-11s %s\n' "$d" "$s" "$(field "$r" command)"
  done
  exit 0
fi

if [ "${1:-}" = "--status" ]; then
  owner=$(lock_owner)
  if [ -z "$owner" ]; then echo "No run in progress."; exit 0; fi
  pid=${owner%% *}; id=${owner#* }
  echo "Running: $id (worker pid $pid; stop with: kill $pid)"
  sed 's/^/  /' "$runs/$id/record.txt"
  echo "--- last lines of $runs/$id/log.txt"
  tail -n 15 "$runs/$id/log.txt" 2>/dev/null
  exit 0
fi

# ----------------------------------------------------------------- launch ---
threads=${JULIA_NUM_THREADS:-$default_threads}; pull=1
while [ $# -gt 0 ]; do
  case "$1" in
    -t) [ $# -ge 2 ] || die "-t needs a value"; threads=$2; shift 2 ;;
    --no-pull) pull=0; shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; break ;;
    -*) die "unknown option $1" ;;
    *) break ;;
  esac
done
[ $# -ge 1 ] || die "usage: run.sh [-t THREADS] [--no-pull] SCRIPT [ARGS...]"
script=$1; shift

nthreads=$(total_threads "$threads") \
  || die "threads must be a number or N,M, not '$threads' ('auto' would use every core; the limit is $max_threads)"
[ "$nthreads" -ge 1 ] || die "threads must be at least 1, got '$threads'"
[ "$nthreads" -le "$max_threads" ] || die "threads '$threads' exceed the limit of $max_threads"
# BLAS threads from the script's own `blas_threads=K` argument (the last one wins; the
# scripts default to 1, which runs BLAS on the calling thread and adds none)
blas=1
for arg in "$@"; do
  case "$arg" in blas_threads=*) blas=${arg#blas_threads=} ;; esac
done
[[ $blas =~ ^[0-9]+$ ]] && [ "$((10#$blas))" -ge 1 ] \
  || die "blas_threads must be a whole number of at least 1, not '$blas'"
blas=$((10#$blas))
[ "$blas" -le 1 ] || [ $((nthreads + blas)) -le "$max_threads" ] \
  || die "-t $threads plus blas_threads=$blas exceeds the limit of $max_threads threads"

cd "$repo" || die "cannot cd to $repo"
mkdir -p "$runs"

# Lock first, so a pull can't change the code under a running job.
if ! mkdir "$lock" 2>/dev/null; then
  owner=$(lock_owner)
  [ -z "$owner" ] || die "busy: run ${owner#* } (worker pid ${owner%% *}) is still running"
  echo "run.sh: removing stale lock" >&2
  rm -rf "$lock"; mkdir "$lock" || die "cannot take lock $lock"
fi
echo $$ >"$lock/pid"   # the launcher holds it until the worker takes over
release() { rm -rf "$lock"; }

git diff --quiet HEAD -- || { release; die "tracked files have local changes; commit or stash them first"; }
if [ "$pull" -eq 1 ]; then
  git pull --ff-only || { release; die "git pull --ff-only failed"; }
fi
[ -f "$script" ] || { release; die "no such script: $script"; }

id="$(date +%Y%m%d-%H%M%S)-$(basename "$script" .jl)"
rundir="$runs/$id"
mkdir "$rundir" || { release; die "run directory $rundir already exists"; }
echo "$id" >"$lock/run"

{
  echo "run: $id"
  echo "command: julia --project=. -t $threads $script${*:+ $*}"
  echo "commit: $(git rev-parse HEAD)"
  echo "branch: $(git rev-parse --abbrev-ref HEAD)"
  echo "pulled: $([ "$pull" -eq 1 ] && echo yes || echo no)"
  echo "host: $(hostname)"
  echo "cores: $(ncores)"
  echo "threads: $threads"
  echo "blas threads: $blas"
  echo "julia: $("${JULIA:-julia}" --version 2>/dev/null || echo unknown)"
  echo "start: $(date '+%Y-%m-%d %H:%M:%S %z')"
  echo "log: output/runs/$id/log.txt"
} >"$rundir/record.txt"

nohup "$self" __worker "$rundir" "$threads" "$script" "$@" </dev/null >/dev/null 2>&1 &
wpid=$!
# Wait for the worker to take the lock over, so it is never held by a dead pid.
for _ in $(seq 50); do
  [ "$(cat "$lock/pid" 2>/dev/null)" = "$wpid" ] || [ ! -d "$lock" ] && break
  sleep 0.1
done

echo "Started $id (worker pid $wpid)."
echo "  log:    output/runs/$id/log.txt"
echo "  status: scripts/remote/run.sh --status"
