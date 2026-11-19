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
