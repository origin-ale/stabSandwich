# Item 8: τ for the `:chi3` layer-end disentangler

Post-port plan item 8, decision 1: the layer-end disentangler drops only numerical zeros, i.e.
singular values σ ≤ τ·σ_max on the bond. These scripts measure where the search candidates'
singular values lie, so that τ can be placed in the gap between rounding noise and physical
values.

- `tau_measure.jl N p1,p2,... samples out.jls`: replays the `:camps` layer loop of
  `Comptime_prob_tm.jl` (same seeds and doping, `thl = 1e-12`) and records σ/σ_max of the
  identity and the 19 candidate gates on every bond of dimension ≠ 1, at each layer end,
  before the disentangler. See its header for what exactly is recorded.
- `tau_analyse.jl out.jls ...`: histogram, widest empty gap, and per-τ counts (see its header).

Run from the repository root, one Julia process at a time:

```
julia -t 1 --project=. scripts/item8/tau_measure.jl 12 0.01,0.02,0.04,0.1 10 $TMPDIR/tau_N12.jls   # ~40 s
julia -t 1 --project=. scripts/item8/tau_measure.jl 46 0.01 3 $TMPDIR/tau_N46.jls                 # ~1 min
julia --project=. scripts/item8/tau_analyse.jl $TMPDIR/tau_N12.jls $TMPDIR/tau_N46.jls
```

## Result (2026-11-19, code before item 8, where `:chi3` never acts)

| | Largest rounding value | Smallest physical value | Empty gap |
|---|---|---|---|
| N=12, p ∈ {0.01, 0.02, 0.04, 0.1}, samples 1–10 (689 searches) | 6.2e-15 | 2.0e-5 | ~9.5 decades |
| N=46, p = 0.01, samples 1–3 (1324 searches) | 5.1e-15 | 2.7e-2 | ~12.7 decades |

- Chosen: **τ = 1e-10**. No candidate value lies within a factor of 100 of it.
- Some gate lowers the τ-rank in 269/689 (N=12) and 239/1324 (N=46) searches, the same for
  every τ from 1e-14 to 1e-7; at τ = 1e-16 the count drops (too tight: zeros are missed).
- The `truncated` flag never flipped, so the disentangler runs at every layer end.
- The stored bonds carry no zero padding (every stored dimension equals the numerical rank at
  τ = 1e-10). The zeros in the identity's spectrum come from the two-site SVD returning
  min(2χ_L, 2χ_R) values.

## After the fix

Item 8's Verification repeats this check on the final code: rerun both commands; the
"vals in [τ/100,100τ]" column must be 0 at τ = 1e-10. `tau_measure.jl` uses CampsPP internals
(`_ID_2Q`, `BondSearchDisentangler(...).matrices`, `layer_disentangler`, `layer_disentangle!`),
so adapt it if item 8 renames them. It records only the first search on each bond at a layer
end, not the searches on a state already changed by a gate earlier in the same sweep.

## Verification on the final code (2026-11-20, `989f6b6`, laptop, 1 thread)

**τ rerun: PASS.** Same two commands as above. At τ = 1e-10 the "vals in [τ/100,100τ]"
column is 0 at both N (it stays 0 from 1e-12 to 1e-8 at N=12, and from 1e-13 to 1e-4 at
N=46).

| | Largest rounding value | Smallest physical value | Empty gap |
|---|---|---|---|
| N=12 (320 searches) | 2.6e-15 | 8.1e-6 | ~9.5 decades |
| N=46 (190 searches) | 8.1e-16 | 6.3e-2 | ~13.9 decades |

- `:chi3` acts. N=12: in 25 of the 32 doped samples, 270 gates in all (the 8 undoped samples
  stay product states). N=46: in all 3 samples (87, 21 and 63 gates).
- The `truncated` flag never flips. No stored bond carries zero padding at the next layer end.
  That is weak evidence for `maxdim = r`, because `apply!` has re-SVD'd the bonds by then. The
  "no `maxdim = r`" mutant (caught by 79 checks) is the real evidence.
- Scope: as before, only the first search per bond at a layer end is recorded. The searches
  later in the same sweep, those in later sweeps, and the values that `maxdim = r` drops are
  not recorded.
- There are fewer recorded searches than before item 8 (320 vs 689, and 190 vs 1324). That
  follows from more bonds being back at dimension 1. It is not the number of searches the
  disentangler performs, which now includes extra sweeps and was not measured.

**Accuracy against dense (N=12).** `accuracy_dense.jl` runs `campssrc_circuit_dynamics` and the
dense reference on the same samples as above. `accuracy_compare.jl` compares two of its outputs.
The reference is the commit before the item 8 code (`3bca35a`), unpacked with `git archive`
and run with `--project=<that tree>`:

```
julia -t 1 --project=.      scripts/item8/accuracy_dense.jl 12 0.01,0.02,0.04,0.1 10 $TMPDIR/acc_new.jls
julia -t 1 --project=<tree> scripts/item8/accuracy_dense.jl 12 0.01,0.02,0.04,0.1 10 $TMPDIR/acc_old.jls
julia --project=. scripts/item8/accuracy_compare.jl $TMPDIR/acc_old.jls $TMPDIR/acc_new.jls
```

- max |⟨O⟩ − dense| over all 40 samples and layers: 2.7e-14 before item 8, 1.3e-13 after.
- The new error is larger in 21/40 samples, by up to 26×. Every one of the 21 is a sample where
  `:chi3` acted, and the error grows layer by layer after the first gates (e.g. p = 0.1, s = 2:
  3.6e-15 → 6.3e-14 → 1.3e-13). It stays at rounding level, about 1e-14 relative to ‖O‖ ≈ 6.
  Why it grows was not tested: extra SVDs, the changed frame and different `apply!` SVDs could
  each explain it.
- **Criterion, as reworded by the user (2026-11-20): |⟨O⟩ − dense| ≤ 1e-12 at every record
  point, for both versions. PASS** (max 1.3e-13 new, 2.7e-14 old). The original criterion
  ("no larger than old, up to the accumulated truncation error") was not met taken literally.
  It made no allowance for rounding, and the truncation error was not measured.
  The `truncated` flag never flips (seen in the `tau_measure.jl` replays of the same samples).
  That only bounds the total discarded weight below 1e-12: `apply!` still truncates at
  `cutoff = thl`.
- The layer-end max bond dimension (after the disentangler) is never above the old run's
  (0/40 samples). The sum of the per-sample peak bond dimensions drops from 128 to 69.
- The old-vs-new difference was not computed directly. It is bounded by old error + new error,
  ≤ 1.55e-13, well within δ = 1e-6.
- Values are recorded before the layer-end disentangler, so the last layer's disentangle is not
  checked against dense.

**Same check with no truncation (`thl = 0`, 2026-11-20).** This rules out truncation as the
cause of the extra error. `noflag` mode is needed here: at `thl = 0` rounding alone trips the
`truncated` flag in every doped sample (`flag = true`), which would switch the disentangler off.

```
julia -t 1 --project=.      scripts/item8/accuracy_dense.jl 12 0.01,0.02,0.04,0.1 10 $TMPDIR/acc_new_thl0.jls 0 noflag
julia -t 1 --project=<tree> scripts/item8/accuracy_dense.jl 12 0.01,0.02,0.04,0.1 10 $TMPDIR/acc_old_thl0.jls 0 noflag
julia --project=. scripts/item8/accuracy_compare.jl $TMPDIR/acc_old_thl0.jls $TMPDIR/acc_new_thl0.jls
```

- Fidelity check: `noflag` at `thl = 1e-12` on HEAD reproduces the `campssrc` run bitwise
  (errors and layer-end bond dimensions in all 40 samples, 270 gates).
- `apply!` drops nothing at `cutoff = 0`. Before item 8, nothing is dropped anywhere. After
  item 8, the only values dropped are the disentangler's σ ≤ τ·σ_max.
- max |⟨O⟩ − dense|: 2.1e-13 before item 8, 1.5e-13 after. The new error is above the old one
  in 15/40 samples (up to 47×) and below it in others (p = 0.1: s = 1, 4, 5, 6).
- Conclusion: the extra error at `thl = 1e-12` is not a truncation effect. With nothing
  truncated, both versions scatter at rounding level, and which one is larger varies from
  sample to sample. A per-sample "no larger than old" comparison of two rounding errors is
  not meaningful here. Whether the scatter comes from the disentangler's SVDs, its dropped
  near-zeros or the changed frame was not separated. All of these act at the rounding level:
  dropping a physical value would show up at ≳ 1e-10 or more.
- Side effect of `cutoff = 0`: ITensors keeps exact zeros, so the bonds fill with zero padding.
  The layer-end max bond dimension reaches 1999 before item 8 and 256 after, while the exact
  maximum at N = 12 is 64. For this reason the bond dimensions of these runs are not a
  production comparison. The new run is above the old one at some layer end in 4/40 samples
  (e.g. p = 0.02, s = 4: peak 5 → 18), and the sum of the peaks drops from 7869 to 501.

**T11 B replay (local, 1 thread): PASS.** Output is in `output/t11/item8-local/B`. 20 cases,
max |Δ| 3.0e-12. That is about 10× the item 4 replay's 3.5e-13 (`output/t11/item4-local/B`):
the same rounding growth as at N=12, still far below δ = 1e-6. Peak bond dimension at p = 0.01 drops: the old values 16, 2, 16, 32, 2, 8, 4,
40, 8, 16 become 2, 1, 2, 2, 1, 8, 2, 18, 4, 3. At p = 0 it stays at 1.
