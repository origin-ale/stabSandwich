# Phase 5: benchmarks and timings

Local micro-benchmarks of the pp → Pauli-sum converter and of `expectation`, old stack vs new,
plus the remote timing rerun at an equal thread count. Plan: Phase 5 in the port plan
("Converter at large sums", "Pauli-sum handling", and what T11 handed on).

| File | Where | What |
|---|---|---|
| `common.jl` | – | helpers shared by the scripts below (T8 state, sum files) |
| `gen_sums.jl` | local | real N = 46 PP sums of 10³–10⁵ terms, written to `output/phase5/sums/` |
| `bench_new.jl` | local | new stack: converter, conjugation, `expectation` seam (MPO) and a per-term prototype |
| `bench_old.jl` | local, old project | old stack: converter (first and second call per size), old `expectation`, term by term |
| `timing_summary.jl` | local | per-p summary of the per-sample times in `compare_prob_tm.jl` tables |
| `benchenv/` | – | environment holding only BenchmarkTools, stacked on the repo's by `bench_new.jl` |

## Local benchmarks

All on one Julia thread and one BLAS thread, under a timeout (`perl -e 'alarm shift; exec @ARGV' <s> …`;
macOS has no `timeout`). Run from the repository root:

```bash
julia -t 1 --project=. scripts/phase5/gen_sums.jl E=1620                 # 1161 terms
for t in 1e-3 1e-4 1e-5 1e-10; do
  julia -t 1 --project=. scripts/phase5/gen_sums.jl E=2331 thl=$t        # 3459, 16318, 82321, 350864 terms
done
julia -t 1 --project=. scripts/phase5/bench_new.jl                       # → output/phase5/bench_new.tsv
julia -t 1 --project=../stabSandwich-old scripts/phase5/bench_old.jl     # → output/phase5/bench_old.tsv
```

**Sums.** All are PP-alone outputs on the T8 sample's circuit (p = 0.06, sample 1). Only
`sum_0201` is a physical hand-off, so on the T8 state the other sums are benchmark inputs, not
physical values. Layer ends give 1161–1395 terms up to gate 2295, then 350,864 from gate 2331
on: the triple ending at 2331 multiplies the sum about 250-fold, and inside it gate 2329 gives
1.6·10⁷ terms. The 10⁴ and 10⁵ sizes therefore come from a larger PP threshold at gate 2331
(thl = 1e-3, 1e-4, 1e-5; the production value 1e-10 gives the 350,864-term sum).

**State.** The T8 hand-off state (switch gate 1676, bond dimension 112), rebuilt on each stack.

## Remote timing rerun (you, on `carlos`)

Old and new at the same thread count (`-t 10`; BLAS gets as many again). The T11 replays at
`-t auto` = 192 can't be compared with anything run under the cap.

Commit and push this directory first, then, on `carlos`:

```bash
cd campspp/stabSandwich && git pull               # port-magictensors
# Old code: a worktree at the pre-port commit (once). Its run.sh predates the thread cap, so
# always pass -t; --no-pull because the worktree is on a detached HEAD.
git worktree add ../stabSandwich-old 4c49e78
cd ../stabSandwich-old && julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 1. Old A, then old B (each its own run; wait for the first to finish). `p=` is the old
#    script's argument; the current Comptime_prob_tm.jl takes `magic_probs=` instead
scripts/remote/run.sh --no-pull -t 10 Comptime_prob_tm.jl methods=campspp p=0,0.01,0.02,0.04 samples=25
scripts/remote/run.sh --no-pull -t 10 Comptime_prob_tm.jl methods=camps p=0,0.01 samples=10

# 2. New A and B (main clone), replaying the Phase 0 inputs
cd ../stabSandwich
scripts/remote/run.sh -t 10 scripts/t11/replay_prob_tm.jl run=output/runs/20261110-150852-Comptime_prob_tm out=output/phase5/t10/A
scripts/remote/run.sh -t 10 scripts/t11/replay_prob_tm.jl run=output/runs/20261110-160640-Comptime_prob_tm out=output/phase5/t10/B

# 3. Thread scan of the B outlier (p = 0.01, sample 8) on the new stack; -t 10 is in run 2
scripts/remote/run.sh -t 1 scripts/t11/replay_prob_tm.jl run=output/runs/20261110-160640-Comptime_prob_tm methods=camps p=0.01 samples=8 out=output/phase5/threads/t1
scripts/remote/run.sh -t 4 scripts/t11/replay_prob_tm.jl run=output/runs/20261110-160640-Comptime_prob_tm methods=camps p=0.01 samples=8 out=output/phase5/threads/t4
```

**One run at a time across both clones.** The worktree has its own `output/runs/.lock`, so
`run.sh` won't stop a run in one clone while the other is busy; check `run.sh --status` in both.

Fetch all six runs, locally. The two old runs (A, B) come from the worktree, so they need
`STAB_DIR`. The four new runs (A, B, scan t1, scan t4) come from the main clone, which is the
default:

```bash
for id in <old A> <old B>; do STAB_DIR=campspp/stabSandwich-old scripts/remote/fetch.sh $id; done
for id in <new A> <new B> <scan t1> <scan t4>; do scripts/remote/fetch.sh $id; done
```

The old script writes fixed names (`output/comptime_prob_tm.txt`, …), which the Phase 0 runs and
then old A already put there locally. So `fetch.sh` may report conflicts and leave the local
copies alone. That's expected and harmless: don't pass `-f`. The comparisons below read each
run's own copy under `output/runs/<id>/files/`.

Compare and summarise, locally (the comparison is also a fresh old-vs-new regression check):

```bash
julia --project=. scripts/t11/compare_prob_tm.jl old=output/runs/<old A>/files/output \
  new=output/runs/<new A>/files/output/phase5/t10/A tsv=output/phase5/timing_A.tsv
julia --project=. scripts/t11/compare_prob_tm.jl old=output/runs/<old B>/files/output \
  new=output/runs/<new B>/files/output/phase5/t10/B tsv=output/phase5/timing_B.tsv
for t in 1 4; do
  julia --project=. scripts/t11/compare_prob_tm.jl old=output/runs/<old B>/files/output \
    new=output/runs/<scan t$t>/files/output/phase5/threads/t$t tsv=output/phase5/timing_B_t$t.tsv
done
julia scripts/phase5/timing_summary.jl output/phase5/timing_*.tsv
```

## Results

### Local benchmarks (2026-11-17, MacBook Air, Julia 1.13.3, 1 thread, BLAS 1 thread)

Raw tables: `output/phase5/bench_new.tsv` and `bench_old.tsv`; logs in `output/phase5/logs/`.
The state is the T8 hand-off state, rebuilt on each stack: old 142 s, new 23 s, both switching at
gate 1676 with k = 37 and bond dimension 112. New-stack times are BenchmarkTools medians; the MPO
seam and per-term runs estimated above 60 s are single runs. Old-stack times are single
`@elapsed` runs.

| terms | old conv. (1st call) | old conv. (2nd) | new conv. | old MPO | new MPO seam | new MPO maxdim | old term by term | new per-term, naive | new per-term, prefix |
|---|---|---|---|---|---|---|---|---|---|
| 201 | 1.53 s | 0.48 ms | 0.07 ms | 1.97 s | 1.63 s | 99 | 5.2 s | 1.37 s | 0.82 s |
| 538 | 5.18 s | 2.7 ms | 0.19 ms | 21.2 s | 20.4 s | 381 | 13.9 s | 3.75 s | 3.48 s |
| 1161 | 12.4 s | 12.7 ms | 0.40 ms | 280 s | 282 s | 904 | 29.9 s | 8.30 s | 7.70 s |
| 3459 | – | – | 1.2 ms | – | – | – | – | 25.9 s | 13.2 s |
| 16,318 | – | – | 5.9 ms | – | – | – | – | 123 s | 45.8 s |
| 82,321 | – | – | 32 ms | – | – | – | – | 626 s | 155 s |
| 350,864 | – | – | 0.22 s | – | – | – | – | (est. 2500 s, not run) | (est. 2700 s from an unsorted subset, not run) |

- **Values:** old MPO vs new seam: 2.4e-13, 7.3e-14 and 6.4e-15 at 201, 538 and 1161 terms; old
  MPO vs new per-term ≤ 2.7e-13; old MPO vs the T8 reference 4.1e-13 and 1.7e-13; the two
  per-term variants agree with each other to ≤ 3e-14 at every size.
- **Old converter: per-size compilation confirmed.** The first call at a new term count costs
  1.5 s (201 terms) to 12.4 s (1161), and a second call 0.5–13 ms. The probe tells new content
  apart from a new size: a *different* 201-term sum converts in 0.58 ms (a 201-term method already
  exists), while 202 and 203 terms take 1.56 s each and a repeat of 202 takes 0.51 ms. So every
  old hand-off or recorded ⟨O⟩ at a term count not seen before in the process paid
  about 1.5–12 s of compilation. The new converter is linear (about 0.6 µs/term) and compiles
  nothing per size.
- **Conjugation:** new 0.13 ms (201 terms) to 0.41 s (350,864 terms), linear and negligible.
- **MPO seam (old and new alike):** most of the time goes into the MPO build and `inner`, and the
  cost grows with the MPO bond dimension (99 → 381 → 904). 538 → 1161 terms costs 14×, about
  terms^3.4 on these sums. Old and new are equally fast (same ITensorMPS route).
- **Per-term prototype:** linear. Prefix sharing (terms sorted by Pauli string, left
  environments reused) costs about 2 ms/term at 82k terms, 4× faster than naive there. It already
  beats the MPO at 201 terms (0.8 s vs 1.6 s), and at 1161 terms it is 37× faster (7.7 s vs 282 s).
  It isn't in `src/`; replacing the MPO route is the optional decision in the plan.

### Remote timing rerun (2026-11-17, `carlos`, Julia 1.13.3, `-t 10`, BLAS threads = Julia threads)

| run | id | wall time |
|---|---|---|
| old A (`4c49e78`) | `20261117-130003-Comptime_prob_tm` | 24.5 min |
| old B (`4c49e78`) | `20261117-132531-Comptime_prob_tm` | 11.3 min |
| new A (`d6bdfe7`) | `20261117-133845-replay_prob_tm` | 3.7 min |
| new B (`d6bdfe7`) | `20261117-134249-replay_prob_tm` | 6.1 min |
| new B outlier, `-t 1` | `20261117-134935-replay_prob_tm` | 2.0 min |
| new B outlier, `-t 4` | `20261117-135245-replay_prob_tm` | 2.1 min |

The six ran one after another, with no overlap. Comparison tables: `output/phase5/timing_*.tsv`;
log: `output/phase5/logs/compare.log`.

**Regression check: PASS, 0 violations** (A 100 cases, B 20, scans 1 + 1). Switch gate, k, peak
bond dimension and peak n. of Paulis are identical on every sample. Max |Δ⟨O⟩| is 8.5e-13 on CAMPS
layers and 1.3e-12 on PP layers (A), and 3.2e-13 on B.

**Times** (`timing_summary.jl`; indicative, shared cluster):

| job | p | mean old | mean new | median × | min × |
|---|---|---|---|---|---|
| A (CAMPS-PP, 25 samples) | 0 | 13.4 s | 1.05 s | 11.8 | 10.5 |
| | 0.01 | 13.9 s | 1.27 s | 11.5 | 8.2 |
| | 0.02 | 14.0 s | 1.44 s | 11.0 | 5.6 |
| | 0.04 | 12.3 s | 1.95 s | 6.5 | 3.7 |
| B (CAMPS, 10 samples) | 0 | 20.6 s | 6.55 s | 2.8 | 2.5 |
| | 0.01 | 33.8 s | 19.6 s | 1.8 | 1.3 |

Every ratio is higher than in the T11 replays at `-t auto` = 192 (A 10/9/8/3×, B 2×/~1×).

**Thread scan, B outlier (p = 0.01, sample 8), new stack:** 29.0 s at `-t 1`, 38.6 s at `-t 4`,
43.9 s at `-t 10`, against 722 s at `-t auto` = 192 (T11) and 19.7 s on the laptop at 1 thread.
The old stack takes 58.5 s at `-t 10`. The time grows with the thread count. This supports the
hypothesis that BLAS/Strided threading oversubscribes the small dense operations; it is one
sample per setting, so it's indicative only. In practice: few threads per process, and
parallelise over samples.
