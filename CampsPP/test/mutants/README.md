# T10 mutants

Mutation tests for the ported (MagicTensors) code paths: 22 hand-picked mutants, each a small
patch to `CampsPP/src`, run against the full fast-tier suite. They are **not** part of
`Pkg.test()`; run them by hand, e.g. after the Phase 6 deletions, after a MagicTensors upgrade,
or after changing the tests they rely on.

- `patches/<ID>.patch`: one mutant each (`patch -p1` relative to the repo root).
- `gen_patches.py`: regenerates the patches from exact-string replacements (the list in the file
  is the definition of each mutant). Use it when a patch no longer applies after `src` changes:
  `python3 -I gen_patches.py <CampsPP dir> patches`.
- `run_mutants.sh`: copies the repo's `CampsPP` to a scratch directory and
  runs one fresh `Pkg.test()` per mutant there (the repo is never modified), with a wall-clock
  cap per run.
- `who_caught.py`: which test files and testsets caught each mutant, and how many of those
  checks are on the old stack or compare against the old code.

```sh
bash CampsPP/test/mutants/run_mutants.sh "$TMPDIR/mutants"            # BASE + all 22, ~1.5-2 h
bash CampsPP/test/mutants/run_mutants.sh "$TMPDIR/mutants" 600 BASE E7 E12
python3 -I CampsPP/test/mutants/who_caught.py "$TMPDIR/mutants/results" E12 T6
```

A run takes about 4 min (fast tier plus precompiling the mutated package). Run one runner at a
time: every run is a full Julia test process.

## The mutants and their 2026-11-16 results

Results on `3384883`, before the hang guard described below. The number is the failing checks.

| ID | File | Mutation | Result |
|---|---|---|---|
| E1 | EvolveQubitCAMPS.jl | `mod(m, 4)` → `mod(m, 2)` (quarter turns) | killed, 93 |
| E3 | EvolveQubitCAMPS.jl | rotate commuting instead of anticommuting generators | killed, 207 |
| E4 | EvolveQubitCAMPS.jl | swap the n = 1 and n = 3 tableaux | killed, 188 |
| E5 | EvolveQubitCAMPS.jl | frame updated with G instead of G† | killed, 188 (≡ E4) |
| E6 | EvolveQubitCAMPS.jl | frame update on the wrong side | killed, 223 |
| E7 | EvolveQubitCAMPS.jl | Heisenberg image through C instead of C† | killed, 259 (hung, see below) |
| E8 | EvolveQubitCAMPS.jl | disentangler written as D instead of D⁻¹ | killed, 167 (hung) |
| E9 | EvolveQubitCAMPS.jl | magic-state angle ignores the sign of Q | killed, 119 |
| E10 | EvolveQubitCAMPS.jl | k never incremented | killed, 364 |
| E11 | EvolveQubitCAMPS.jl | `:logical` branch rotates P instead of Q | killed, 176 |
| E12 | EvolveQubitCAMPS.jl | truncating `apply!` ignores `thl` | killed, 6 (see below) |
| E13 | EvolveQubitCAMPS.jl | non-unitary magic-state matrix | killed, 279 |
| E14 | EvolveQubitCAMPS.jl | rotation exp(−iϕP) instead of exp(iϕP) | killed, 238 |
| D1 | Disentangler.jl | qubit k+1 classified as logical | killed, 599 |
| D3 | Disentangler.jl | sign of Q dropped | killed, 388 |
| D5 | Disentangler.jl | CNOT control and target swapped in `build_D` | killed, 462 (hung) |
| P1 | MagicTensorsPauliPropagation.jl | complex coefficients conjugated | killed, 23 |
| Q2 | DisentangleQubitCAMPS.jl | entropy score uses σ instead of σ² | killed, 302 |
| Q4 | DisentangleQubitCAMPS.jl | search without `orthogonalize!` | killed, 74 (mostly errors, see below) |
| T1 | TimeEvolution.jl | CAMPS→PP switch one gate late | killed, 18 |
| T2 | TimeEvolution.jl | PP layer index off by one after the hand-off | killed, 14 |
| T6 | TimeEvolution.jl | MPS not renormalised | killed, 1 |

Every killed mutant has failing checks outside the old-stack and old-vs-new tests, so the Phase 6
deletions should leave them killed; rerun after Phase 6 to confirm (`who_caught.py` flags any
mutant whose failures are all old).

## Rerun after Phase 6 (2026-11-17)

On `9e5b95c` with the Phase 6 test changes (old stack and old-vs-new tests deleted), cap 600 s:
BASE 5123 pass, 0 fail, 0 broken (145 s). **22 killed, D6 survived (equivalent, as before), no
timeouts or hangs**; every run 135-151 s. Failing checks: D1 413, D3 336, D5 383, E1 107, E3 196
(13 errors), E4/E5 168, E6 203, E7 238, E8 151, E9 107, E10 348, E11 167, E12 6, E13 264, E14 202,
P1 17, Q2 8, Q4 43 (40 errors), T1 9, T2 7, T6 1.

- **Q2 is now thin (302 → 8).** Most of its catches were the deleted old-vs-new disentangler
  parity tests. Left: the `max_iter = N` sweep test (`TestTimeEvolution.jl`) and the Bell-pair
  control where `:entangle` must gain exactly 1 bit (`TestDisentangleQubitCAMPS.jl`). The gate
  choices of `:entangle` are now pinned only through these.
- **Q4 (43, of which 40 errors)** is still caught first by MagicTensors' orthogonality-centre
  check in `apply_and_svd`; the value-level catches are 3 (the frame-change counts in
  `TestTimeEvolution.jl` and the snake skip rule).

## Notes

- **D6 deleted (2026-11-19).** It made `reducetoX` return S†·P·S (+X) instead of S·P·S† (−X).
  It was equivalent (the sign of the reduced string is never read), and post-port item 5 made
  that change in `src`, so the mutant would now be identical to the code. The 2026-11-16 and
  2026-11-17 results above predate the deletion.
- **E4 ≡ E5**: swapping the n = 1 and n = 3 tableaux is the same as applying G instead of G†.
- **E12 is caught only by the two threshold ladders** (`TestApplyDense.jl`, `TestTimeEvolution.jl`),
  since every other test runs at thl ≤ 1e-11. Known limitation (user, 2026-11-16): nothing ties
  the error to the *value* of `thl`, so a mis-scaled threshold (e.g. `thl/10`) would most likely
  survive.
- **Q4 is caught first by MagicTensors**: `apply_and_svd` requires the orthogonality centre on
  the bond. If that check goes upstream, the remaining catch is weaker.
- **Hangs.** E7, E8 and D5 made the bond dimension blow up in the N = 46 hand-off rebuild of
  `TestMagicTensorsPP.jl`, which then used a plain `evolve` with no bond-dimension stop, so the
  suite hung instead of failing. The rebuild now uses `evolve_bonddim` with the reference's own
  χ = 64 stop and asserts the switch gate, which turns those hangs into failures. The runner's
  cap is the second line of defence. Checked 2026-11-16 with this runner: BASE 7561 pass,
  0 fail; E7 killed (262 checks, 265 s), E8 killed (170, 251 s), D5 killed (465, 231 s), no
  hangs.
- **Runner history.** The plan asked for one persistent session (Revise). That runner misfiled
  results (`include` wraps the test exception in a `LoadError`), its liveness check missed
  keyword-argument methods and stale world ages, and Revise's file polling raced with the patches.
  Fresh processes cost ~2 min more per mutant and avoid all of it.
