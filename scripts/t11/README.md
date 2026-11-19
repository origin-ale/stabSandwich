# T11: old-vs-new paired regression

These scripts replay the Phase 0 `Comptime_prob_tm.jl` samples on the MagicTensors
`QubitCAMPS` stack and compare them with the old runs, sample by sample. The plan
is in "Old-vs-new regression design" (Phase T) and "Phase 0" of the port plan.

| File | Where | What |
|---|---|---|
| `replay_prob_tm.jl` | remote (through `run.sh`) or local | rebuilds each sample from its saved bitstring and doped gates, then runs the new stack with the run's own parameters |
| `compare_prob_tm.jl` | local | per-sample comparison and the old/new resource table; exit 1 on any violation |
| `sensitivity.sh` | local | injects a known defect into a few replays and checks that the comparison flags it |

## Phase 0 reference runs

- **A** `output/runs/20261110-150852-Comptime_prob_tm`: `:campspp`, p ∈ {0, 0.01, 0.02, 0.04}, 25 samples
- **B** `output/runs/20261110-160640-Comptime_prob_tm`: `:camps` (`:chi3`/`:snake`), p ∈ {0, 0.01}, 10 samples

## Remote runs (you, on `carlos`)

Replays run from `main`. Commit and push first: `run.sh` refuses a dirty tree, and it
pulls (`git pull --ff-only`) before each run. It doesn't instantiate, so run
`Pkg.instantiate()` on `carlos` after any commit that changes a Manifest (MagicTensors
is fetched from GitHub). The port branch `port-magictensors` was merged and deleted; a
clone still on it can't pull, so switch it to `main` (`git fetch --prune && git checkout
main`).

```bash
cd campspp/stabSandwich
git pull --ff-only
julia --project=. -e 'using Pkg; Pkg.instantiate()'
scripts/remote/run.sh -t 10 scripts/t11/replay_prob_tm.jl run=output/runs/20261110-150852-Comptime_prob_tm out=output/t11/<label>/A
# once A has finished (one run at a time):
scripts/remote/run.sh -t 10 scripts/t11/replay_prob_tm.jl run=output/runs/20261110-160640-Comptime_prob_tm out=output/t11/<label>/B
```

`out=` keeps a new replay apart from earlier ones (the default is
`output/t11/<Phase 0 id>/`); `<label>` names the change being checked, e.g. `item3`.
B runs `:camps` only: PP builds its gates (`pp.PauliRotation` → `getpauli`) but propagates
nothing, so A carries most of the weight for a change to PauliPropagation.

**Threads:** at most 20 at once on `carlos`. `run.sh` refuses `-t auto` and anything
above 20, and defaults to 10. The scripts also set BLAS to the same number of
threads, so a run can use about twice its `-t`: keep `-t` at 10 or below. The Phase 0 runs and the first replays (2026-11-16) used
`-t auto`; times from a different thread setting aren't comparable with those.
Times are indicative only, because the cluster is shared.
Each replay writes to its `out=` directory (default `output/t11/<Phase 0 id>/`), and
`run.sh` copies that into its own run directory. Its log (`log.txt` in the run directory)
records the loaded versions of PauliPropagation, MagicTensors, QuantumClifford, ITensors
and ITensorMPS. Each version comes with the install directory's slug, which identifies
the exact code even when the version number doesn't change, as with MagicTensors.

## Comparison (local, after `fetch.sh`)

```bash
julia --project=. scripts/t11/compare_prob_tm.jl \
  old=output/runs/20261110-150852-Comptime_prob_tm/files/output \
  new=output/runs/<replay A id>/files/output/t11/<label>/A
```

Without `out=`, the replay is under `files/output/t11/20261110-150852-Comptime_prob_tm`.

Do the same for B. Always read the reference from the run directories, never from
`output/`: the Phase 0 jobs share file names there.

Tolerances (user, 2026-11-10):
- CAMPS-PP layers that both stacks ran in CAMPS: 1e-10.
- Everything else: δ = 1e-6.

The CAMPS-PP switch gate, k at the switch and the PP stop gate must be identical,
and so must the inputs (bitstrings and doped gates). A sample above tolerance is
investigated individually, never re-rolled or averaged away. δ is never loosened
without your written approval. Bond dimension, n. of Paulis and times are printed
for review and are not pass/fail.

## Optional: `Randrot_complexity.jl` (local, N = 12)

`randrot_complexity.jl` runs the script's computation without the plotting, on
either stack. The old stack runs from the pre-port worktree, and the new stack
from `main`:

```bash
julia --project=../stabSandwich-old scripts/t11/randrot_complexity.jl stack=old out="$TMPDIR/rr_old.dat"
julia --project=. scripts/t11/randrot_complexity.jl stack=new out="$TMPDIR/rr_new.dat"
julia scripts/t11/compare_randrot.jl "$TMPDIR/rr_old.dat" "$TMPDIR/rr_new.dat"
```

M₂/N must agree within |z| ≤ 3 at every point. EE/N depends on the frame, so it is
compared side by side only.

## Sensitivity check (local)

```bash
scripts/t11/sensitivity.sh output/runs/20261110-150852-Comptime_prob_tm "$TMPDIR/t11"
```

The check runs a clean replay (which must pass) and two defective ones: the magic
phase × 1.05, and every doped pair moved one XXZ triple later. Each defective
replay must be flagged by the ⟨O⟩/switch checks, not only by the input check.
