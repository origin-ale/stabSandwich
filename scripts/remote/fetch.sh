#!/bin/bash
# Fetch a finished remote run into the local repo, verifying checksums.
# Run this LOCALLY, from anywhere inside the repo.
#
#   scripts/remote/fetch.sh RUN-ID     # fetch one run ("latest" = newest finished run)
#   scripts/remote/fetch.sh -f RUN-ID  # also overwrite local files that differ
#   scripts/remote/fetch.sh --list     # list remote runs and their status
#   scripts/remote/fetch.sh --status   # show the remote run in progress
#
# Remote location (environment variables):
#   STAB_HOST  ssh host alias of the target  (default: stab, see ssh_config.example)
#   STAB_DIR   repo clone on the target, relative to the remote home or absolute
#              (default: stabSandwich)
#
# Steps:
#   1. refuses runs that haven't finished (no SHA256SUMS on the target yet),
#   2. streams the run directory (record, log, and the run's own copy of its
#      output files) through one ssh connection (tar) into a staging directory,
#   3. checks every file against SHA256SUMS; on any mismatch nothing is installed,
#   4. installs the run directory as output/runs/RUN-ID/ and copies each output
#      file to its own path under output/, where scripts and plots expect it.
#      A local file there with different content is left alone and reported,
#      unless -f is given. The run directory keeps its copy either way.
#
# A non-zero exit status means something was not installed.
set -uo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd -P)"
host=${STAB_HOST:-stab}
rdir=${STAB_DIR:-stabSandwich}
ssh_cmd=${STAB_SSH:-ssh}   # overridable for testing only

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; }
die() { echo "fetch.sh: $*" >&2; exit 1; }

sha256() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$@"; else shasum -a 256 "$@"; fi
}

q() { printf '%q' "$1"; }   # quote for the remote shell
remote() { "$ssh_cmd" "$host" "cd $(q "$rdir") && $1"; }

force=0
case "${1:-}" in
  --list)   remote "scripts/remote/run.sh --list"; exit $? ;;
  --status) remote "scripts/remote/run.sh --status"; exit $? ;;
  -h|--help|"") usage; [ -n "${1:-}" ]; exit $? ;;
  -f) force=1; shift ;;
esac
[ $# -eq 1 ] || die "usage: fetch.sh [-f] RUN-ID | --list | --status"
id=$1

if [ "$id" = "latest" ]; then
  id=$(remote "ls -1 output/runs 2>/dev/null | LC_ALL=C sort -r | while read d; do [ -f output/runs/\$d/SHA256SUMS ] && { echo \$d; break; }; done; true") \
    || die "cannot reach $host"
  [ -n "$id" ] || die "no finished run on $host"
  echo "Latest finished run: $id"
fi
case "$id" in */*|.*|"") die "bad run id: $id" ;; esac

rel="output/runs/$id"
remote "test -d $(q "$rel")" || die "no run $id on $host (or host unreachable)"
remote "test -f $(q "$rel/SHA256SUMS")" \
  || die "run $id has not finished (no SHA256SUMS yet); see: fetch.sh --status"

cd "$repo" || die "cannot cd to $repo"
mkdir -p output/runs
stage=$(mktemp -d "$repo/output/.fetch-$id.XXXXXX") || die "cannot create staging directory"
trap 'rm -rf "$stage"' EXIT

echo "Transferring $id from $host ..."
remote "tar cf - $(q "$rel")" \
  | tar xf - -C "$stage" \
  || die "transfer failed; nothing installed"

# Verify. The checksum list itself came over the same connection, so also make
# sure it covers every file the run reported.
[ -f "$stage/$rel/SHA256SUMS" ] || die "SHA256SUMS missing from the transfer; nothing installed"
missing=$(cd "$stage/$rel" && while IFS= read -r f; do
  grep -qxF "$(sha256 "files/$f" 2>/dev/null)" SHA256SUMS || echo "$f"; done <files.txt)
[ -z "$missing" ] || die "output files missing or not covered by SHA256SUMS: $missing; nothing installed"
if ! (cd "$stage/$rel" && sha256 -c --quiet SHA256SUMS); then
  die "checksum mismatch (truncated or altered transfer); nothing installed. Re-run fetch.sh."
fi
n=$(wc -l <"$stage/$rel/files.txt" | tr -d ' ')
echo "Verified $n output file(s) plus log and record."

# Install.
if [ -e "$rel" ]; then
  if diff -rq "$stage/$rel" "$rel" >/dev/null 2>&1; then :; else
    [ "$force" -eq 1 ] || die "$rel already exists locally with different content; use -f to replace it"
    rm -rf "$rel"
  fi
fi
[ -e "$rel" ] || mv "$stage/$rel" "$rel"

status=0
while IFS= read -r f; do
  src="$rel/files/$f"
  if [ -e "$f" ] && ! cmp -s "$src" "$f"; then
    if [ "$force" -eq 0 ]; then
      echo "  kept local (differs): $f" >&2; status=1; continue
    fi
    echo "  replaced: $f"
  elif [ -e "$f" ]; then
    echo "  unchanged: $f"; continue
  else
    echo "  new: $f"
  fi
  mkdir -p "$(dirname "$f")" && cp -p "$src" "$f" || { echo "  copy failed: $f" >&2; status=1; }
done <"$rel/files.txt"

echo
sed -nE 's/^(status|exit code|commit|warning): /  \1: /p' "$rel/record.txt"
[ "$status" -eq 0 ] || echo "Some local files differ and were kept; re-run with -f to replace them." >&2
exit $status
