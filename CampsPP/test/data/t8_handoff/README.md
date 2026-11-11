# T8 reference data: old-code ⟨O⟩ on real N = 46 Pauli sums

Saved during T-pre (2026-11-11) with the **old** code (vendored CliffordMPS), for T8:
the new pp → `QubitPauliSum` converter and the `expectation` seam must reproduce these
values to about 1e-10 (plan, "Converter at large sums").

## Files

- `state.txt`: a real CAMPS-PP hand-off state. Sample 1 of the `Comptime_prob_tm.jl`
  ensemble at p = 0.06 (rng_seed = 100, μ = 0.6, `:xy` doping at 3π/16), evolved with CAMPS
  (untruncated) up to its switch gate s = 1676 (χ = 64 reached; bond dimension 112, k = 37
  in the old frame). The file lists the excited qubits and the doped gate indices, so the new
  code can rebuild the same physical state by applying gates 1…s of
  `xxz_circuit(π/4, π/4, 23, 46)` with those gates' angles replaced by 3π/16. ⟨O⟩ is
  frame-invariant, so the new frame does not matter.
- `sum_0201.txt`: the **real hand-off** sum: the transferred magnetization propagated by PP
  (thl = 1e-10) through gates s+1…2464, i.e. what CAMPS-PP would hand to the converter when
  recording after gate 2464. ⟨O⟩ is the physical value at that gate.
- `sum_0538.txt`: a **larger, non-physical** sum: PP alone from t = 0 through gates 1…1504.
  Gate 1504 falls inside an XX-YY-ZZ triple (after the XX gate); it is still a genuine PP
  output with realistic term structure, but its ⟨O⟩ on this state has no physical meaning.
  It is only a converter/seam check.
- `generate.jl`: the script that produced them (run from the repo root with
  `julia --project=.`; see its header).

Sum format: header lines start with `#` (end gate `e`, `nterms`, the old ⟨O⟩ through the
MPO route, the old ⟨O⟩ term by term, the old MPO time); then one row per term:
46-character Pauli string (character j = qubit j, CONVENTIONS.md §1) and the `Float64`
coefficient (`repr`, round-trips exactly). Rows are sorted by string. Reading a file back into
a `pp.PauliSum` reproduces the sum exactly (checked when saving).

## Old values

| file | terms | old ⟨O⟩ (MPO) | MPO vs term-by-term | old MPO time (local, one run) |
|---|---|---|---|---|
| `sum_0201.txt` | 201 | 20.45206878679019 | 5.7e-14 | 3.8 s |
| `sum_0538.txt` | 538 | 9.298301330177678 | 1.2e-14 | 21.6 s |

Imaginary parts are ≤ 5e-15. Times include compilation for that sum size and were measured
on a laptop with other work running; they are indicative only (Phase 5 does the timings).

## Scope decisions (2026-11-11)

- The plan asked for sums of about 10², 5·10² and 10³ terms. The user chose to stop at two
  sums: the converter works term by term, so a 538-term sum exercises the same code paths
  (`UInt96` keys, qubits above 32, signs, coefficients, frame) as a 1000-term one. The
  old-reference ceiling is therefore 538 terms; a 10³-term old timing, if wanted for
  Phase 5, comes from the pre-port worktree.
- If a larger sum is wanted later: in this sample, PP alone gives 1161 terms at the end of
  layer 12 (e.g. gate 1521 or 1620, at the end of a full triple), and the old MPO should take
  a few minutes for it. Use `E=<gate>`; don't let the script bisect inside a layer: the term
  count is not monotone in e within a triple (it reaches 60,201 after one XX gate here).
