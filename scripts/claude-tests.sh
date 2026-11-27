#!/bin/bash
# Launch a Claude Code session for test work (Phase T), with the src/ lockdown
# from .claude/lockdown-tests.json layered on top of the managed sandbox.
# Run this from a terminal; extra arguments are passed through to claude.
set -euo pipefail

repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$repo"

settings=".claude/lockdown-tests.json"
[ -f "$settings" ] || { echo "Missing $repo/$settings" >&2; exit 1; }

# Newest bundled binary from the VS Code extension, else claude on PATH.
bin=$(ls -d "$HOME"/.vscode/extensions/anthropic.claude-code-*/resources/native-binary/claude 2>/dev/null \
      | sort -V | tail -1 || true)
[ -n "$bin" ] || bin=$(command -v claude || true)
[ -n "$bin" ] || { echo "No claude binary found (VS Code extension or PATH)" >&2; exit 1; }

# Tells the staged-plan skill (~/.claude/skills/staged-plan) this is the test session.
export STAGED_PLAN_ROLE=test
echo "Test session: $bin --settings $settings" >&2
exec "$bin" --settings "$settings" "$@"
