# Layer-end disentangler for MagicTensors `QubitCAMPS` (CampsPP/src/DisentangleQubitCAMPS.jl,
# plan Phase 3). The specification here is the old CliffordMPS `disentangle!`: the port must
# make the same choices (plan, "Phase 3 status"). Expected values come from hand-checked tableau
# pins and from frame-invariant physics; never from the new implementation. Until Phase 6
# (2026-11-17) they also came from the old code itself: old-vs-new parity of sweeps, gains,
# schedules and frames on generic states (200/200 identical in the Phase 3 probe). Those tests
# went with the old code; `pre-port:CampsPP/test/DisentangleCAMPS/TestDisentangleQubitCAMPS.jl`
# has them.
#
# Decisions this file encodes:
# - Post-port item 8 (user, 2026-11-18/19; CONVENTIONS.md §2, "Numerical zeros"): a singular
#   value σ ≤ τ·σ_max, τ = 1e-10, counts as zero, in the bond search and in the applied gate.
#   `:chi3` scores (numerical rank)³ and now acts; an applied `:chi3` gate strictly lowers the
#   bond dimension; no sweep leaves a numerical zero on a bond. The old no-op `:chi3` (every
#   value counted, zeros kept) is gone, with no option to keep it (plan item 8, decision 2);
#   its testset was replaced, not kept alongside. Expected values below come from the spec,
#   from hand derivations re-checked here on the pinned tableaux, and from frame-invariant
#   physics.
# - `:entangle` keeps the old tie rule (no ε). On Clifford-structured states the chosen gate
#   depends on rounding, so those states are checked through frame-invariant quantities.
#
# The groups marked [item 8] were written before the code and gated on the hook
# `DisentangleCAMPS.NUMERICAL_ZERO_TAU` (plan "How tests are written"). Item 8 has landed, so
# the gate is gone (item 8 review, 2026-11-20: a renamed constant would otherwise have turned
# every group into `@test_broken` with the suite still passing). They now run
# unconditionally, and "[item 8] τ = 1e-10" fails if the hook is missing.

using ITensors, ITensorMPS
using Random
using LinearAlgebra: diag, norm
import CampsPP
import MagicTensors as mt

const DQ = CampsPP.DisentangleCAMPS

"τ of CONVENTIONS.md §2, from the spec (not read from the code)."
const TAU = 1e-10

# -- Helpers -----------------------------------------------------------------------------

"The gates of a `BondSearchDisentangler` as `CliffordOperator`s, in its search order."
search_tableaux(de::DQ.BondSearchDisentangler) = [CliffordOperator(g) for g in de.gates]

"The bond search inside a `sweep_disentangler(…)`."
bond_search(de::DQ.MinDiffDisentangler) = de.disentangler.local_disentangler

"Every N-qubit Pauli string with sign +1."
function all_pauli_ops(N)
  ops = PauliOperator[]
  for code in 0:(4^N - 1)
    xs = falses(N); zs = falses(N)
    for j in 1:N
      d = (code >> (2 * (j - 1))) & 0x3        # I, X, Y, Z = 0, 1, 2, 3
      xs[j] = d == 1 || d == 2
      zs[j] = d == 2 || d == 3
    end
    push!(ops, PauliOperator(0x0, xs, zs))
  end
  return ops
end

"⟨ψ|P|ψ⟩ for every P in `ops`, on the physical state C|mps⟩ (frame-invariant)."
physical_expectations(ψ::mt.QubitCAMPS, ops) = [mt.expectation(ψ, mt.QubitPauli(P)) for P in ops]

# `frame_dagger(ψ)` (the stored frame C† as a copy) comes from helpers.jl.

"A generic state: a random MPS of bond dimension χ behind a random frame C."
function random_state(rng, N, χ)
  s = siteinds("Qubit", N)
  mps = random_mps(rng, ComplexF64, s; linkdims = χ)
  C = random_clifford(rng, N)
  return mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(C), false)
end

"N-qubit |0…0⟩ with a Bell pair (|00⟩ + |11⟩)/√2 on qubits (b, b+1): bond b has dimension 2."
function bell_mps(s, b)
  ψ = MPS(ComplexF64, s, "0")
  mt.apply!(ψ, Matrix(mt.QubitCliffordGate(tCNOT * (tHadamard ⊗ tId1))), [b, b + 1])
  return ψ
end

"Do two MPS hold the same tensors (same indices, same entries)?"
same_tensors(a::MPS, b::MPS) = length(a) == length(b) &&
  all(inds(a[j]) == inds(b[j]) && norm(a[j] - b[j]) == 0 for j in eachindex(a))

# -- Helpers for item 8 (numerical zeros) ------------------------------------------------

"""
Singular values across bond b of `mps`, largest first: the SVD of the two-site block (b, b+1)
of a copy orthogonalized at b, untruncated. A stored bond of dimension χ_b gives at most χ_b
nonzero values; the rest of the min(2χ_{b-1}, 2χ_{b+1}) values are exact or rounding zeros.
"""
function bond_spectrum(mps::MPS, b)
  m = orthogonalize(mps, b)
  left = b > 1 ? (linkind(m, b - 1), siteind(m, b)) : (siteind(m, b),)
  _, S, _ = svd(m[b] * m[b + 1], left...)
  return sort(diag(S); rev = true)
end

"Numerical rank (CONVENTIONS.md §2): the number of values above τ·σ_max."
numrank(σ) = count(>(TAU * maximum(σ)), σ)

"""
Bonds of `mps` that hold a numerical zero: stored dimension above the numerical rank. The
spec's "no bond holds a singular value ≤ τ·σ_max" is `isempty(padded_bonds(mps))`.
"""
padded_bonds(mps::MPS) = [b for b in 1:length(mps)-1 if linkdims(mps)[b] > numrank(bond_spectrum(mps, b))]

"""
A two-qubit state given by its stabilizer generators (s1, s2): is it a product state after
pinned gate k? Hand rule: G|ψ⟩ is stabilized by the images G s G†, and a two-qubit
stabilizer state is a product iff its stabilizer group {s1', s2', s1's2'} holds a weight-1
string. Uses only the pinned tableaux (`cliffconj`, helpers.jl), not the code.
"""
function pinned_disentangles(k, s1, s2)
  a = cliffconj(pinned_tableau(k), s1)
  b = cliffconj(pinned_tableau(k), s2)
  weight(P) = count(j -> xbit(P)[j] || zbit(P)[j], 1:2)
  return any(==(1), weight.((a, b, a * b)))
end

"""
Is the two-qubit Clifford C a local Clifford, or a local Clifford times SWAP? Then the images
of X₁ and Z₁ act on one and the same qubit, and those of X₂ and Z₂ on one qubit.
"""
function local_or_swaplocal(C)
  supp(P) = [j for j in 1:2 if xbit(P)[j] || zbit(P)[j]]
  s = [supp(cliffconj(C, P)) for P in (P"XI", P"ZI", P"IX", P"IZ")]
  return length(s[1]) == 1 && s[1] == s[2] && length(s[3]) == 1 && s[3] == s[4]
end

"""
The pinned gates j that undo pinned gate k on every product state: G_j G_k is local, or local
times SWAP (`A * B` applies B first, TestGoldenSet.jl). On a generic product state |ab⟩
(random single-qubit states), G_j G_k|ab⟩ is a product exactly for these j (with
probability 1), so `:chi3` applies the first of them to G_k|ab⟩.
"""
undoers(k) = [j for j in 1:19 if local_or_swaplocal(pinned_tableau(j) * pinned_tableau(k))]

"N-qubit MPS with amplitudes ψ (DenseReference layout: qubit 1 is the least significant bit)."
mps_from_amplitudes(s, ψ) = MPS(ITensor(reshape(ComplexF64.(ψ), ntuple(_ -> 2, length(s))), s...), s)

"""
N = 4: Φ⁺ = (|00⟩ + |11⟩)/√2 on qubits (2, 3) times (|00⟩ + ε|11⟩)/√(1 + ε²) on qubits
(1, 4). Bond 2 cuts both pairs, so its Schmidt values are {1, 1, ε, ε}/√(2(1 + ε²)): ratio ε
to σ_max. A gate on (2, 3) acts on the Φ⁺ factor only: across bond 2 the spectrum becomes
spec(G Φ⁺) ⊗ {1, ε}/√(1 + ε²).
"""
function bell_times_weak_pair(s, ε)
  ψ = zeros(ComplexF64, 16)
  for b23 in (0, 1), b14 in (0, 1)   # Φ⁺: b2 = b3; weak pair: b1 = b4
    ψ[1 + b14 + 2b23 + 4b23 + 8b14] = (b14 == 0 ? 1 : ε) / sqrt(2 * (1 + ε^2))
  end
  return mps_from_amplitudes(s, ψ)
end

"""
N = 2m: Φ⁺ on qubits (m, m+1) times χ = Σᵢ μᵢ|i⟩|i⟩ (normalised) on the other qubits, with
|i⟩ the basis state of bits i on qubits 1…m−1 and on qubits m+2…N (bit j on qubit j and on
m+1+j). Bond m cuts both factors: its Schmidt values are {1, 1} ⊗ μ, up to the norm, and a
gate on (m, m+1) that makes Φ⁺ a product leaves μ. At N = 4, μ = [1, ε] is
`bell_times_weak_pair`.
"""
function bell_times_schmidt(s, μ)
  N = length(s); m = N ÷ 2
  @assert iseven(N) && length(μ) <= 2^(m - 1)
  ψ = zeros(ComplexF64, 2^N)
  for (i, μi) in enumerate(μ), c in (0, 1)
    ψ[1 + (i - 1) + c * (2^(m - 1) + 2^m) + (i - 1) * 2^(m + 1)] = μi
  end
  return mps_from_amplitudes(s, ψ ./ norm(ψ))
end

"""
Ranks under rules that the spec rejects, for sensitivity checks: `rank_absolute` (σ > τ),
`rank_vs_norm` (σ > τ·√Σσ²) and `rank_cumulative(σ, ref)` (drop the smallest values while
their Σσ² ≤ τ²·ref, as ITensors' `cutoff = τ²` does with ref = Σσ²).
"""
rank_absolute(σ) = count(>(TAU), σ)
rank_vs_norm(σ) = count(>(TAU * norm(σ)), σ)
function rank_cumulative(σ, ref)
  w = cumsum(abs2.(sort(σ)))                    # discarded weight, smallest values first
  return length(σ) - count(<=(TAU^2 * ref), w)
end

"""
A local bond search that logs, for every gate it applies, the bond, its stored dimension
before and after, and its numerical rank after. Wraps the search of `sweep_disentangler`, so
the properties are checked at every applied gate inside real sweeps.
"""
struct RecordingSearch{L <: mt.AbstractLocalDisentangler} <: mt.AbstractLocalDisentangler
  inner::L
  log::Vector{NamedTuple{(:bond, :before, :after, :rank_after), NTuple{4, Int}}}
end
RecordingSearch(inner) =
  RecordingSearch(inner, NamedTuple{(:bond, :before, :after, :rank_after), NTuple{4, Int}}[])

function mt.disentangle!(r::RecordingSearch, ψ::mt.QubitCAMPS, n::Int, right::Bool)
  before = linkdims(mt.get_mps(ψ))[n]
  res = mt.disentangle!(r.inner, ψ, n, right)
  if res[1] > 0
    mps = mt.get_mps(ψ)
    push!(r.log, (bond = n, before, after = linkdims(mps)[n], rank_after = numrank(bond_spectrum(mps, n))))
  end
  return res
end

"`sweep_disentangler(crit, strat, max_iter; min_diff)` with its bond search logged; returns (de, log)."
function recording_sweep(crit, strat, max_iter; min_diff = 0.0)
  rec = RecordingSearch(bond_search(DQ.sweep_disentangler(crit, strat, max_iter)))
  sweep = strat == :snake ? DQ.SnakeDisentangler(rec) : DQ.BrickworkDisentangler(rec)
  return DQ.MinDiffDisentangler(sweep, max_iter; min_diff), rec.log
end

"""
Random product state (random single-qubit states) with pinned gate ks[i] applied on bond
bonds[i] (disjoint qubit pairs), behind a random frame.
"""
function structured_state(rng, N, bonds, ks)
  s = siteinds("Qubit", N)
  mps = random_mps(rng, ComplexF64, s; linkdims = 1)
  for (b, k) in zip(bonds, ks)
    mt.apply!(mps, Matrix(mt.QubitCliffordGate(pinned_tableau(k))), [b, b + 1])
  end
  return mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(random_clifford(rng, N)), false)
end

"A random nonempty set of bonds in 1:N-1, no two adjacent (disjoint qubit pairs)."
function disjoint_bonds(rng, N)
  while true
    bonds = Int[]
    b = rand(rng, 1:2)
    while b <= N - 1
      if rand(rng) < 0.6
        push!(bonds, b); b += 2
      else
        b += 1
      end
    end
    isempty(bonds) || return bonds
  end
end

"Every N-qubit Pauli string up to N = 5; beyond, `nsample` random ones."
function test_paulis(rng, N; nsample = 300)
  N <= 5 && return all_pauli_ops(N)
  return [PauliOperator(0x0, rand(rng, Bool, N), rand(rng, Bool, N)) for _ in 1:nsample]
end

# -- Gate set ----------------------------------------------------------------------------

# The 19 tableaux of the old `Clifford2IndexSet(:entangle)`, in the old order, as the images
# of (X₁, X₂, Z₁, Z₂). Generated from the old `cmps.clifford_gate` on 2026-11-13 and kept
# here so the test survives the removal of `dmrgCliffordMPS`. Hand checks: gate 1 is SWAP
# (X₁ → X₂, …); gate 11 is CNOT with control 1 (X₁ → X₁X₂, X₂ → X₂, Z₁ → Z₁, Z₂ → Z₁Z₂);
# gates 2–10 are iSWAP·(A⊗A) and 11–19 CNOT·(A⊗A) with A ∈ {I, H, H·S}, as the old
# `entangle_set` numbers say.
const OLD_ENTANGLE_PINS = [
  (P"IX", P"XI", P"IZ", P"ZI"),
  (P"ZY", P"YZ", P"IZ", P"ZI"),
  (P"ZY", P"ZI", P"IZ", P"YZ"),
  (P"ZY", P"XZ", P"IZ", P"YZ"),
  (P"IZ", P"YZ", P"ZY", P"ZI"),
  (P"IZ", P"ZI", P"ZY", P"YZ"),
  (P"IZ", P"XZ", P"ZY", P"YZ"),
  (P"ZX", P"YZ", P"ZY", P"ZI"),
  (P"ZX", P"ZI", P"ZY", P"YZ"),
  (P"ZX", P"XZ", P"ZY", P"YZ"),
  (P"XX", P"IX", P"ZI", P"ZZ"),
  (P"XX", P"ZZ", P"ZI", P"IX"),
  (P"XX", P"ZY", P"ZI", P"ZZ"),
  (P"ZI", P"IX", P"XX", P"ZZ"),
  (P"ZI", P"ZZ", P"XX", P"IX"),
  (P"ZI", P"ZY", P"XX", P"ZZ"),
  (P"-YX", P"IX", P"XX", P"ZZ"),
  (P"-YX", P"ZZ", P"XX", P"IX"),
  (P"-YX", P"ZY", P"XX", P"ZZ"),
]

pinned_tableau(i) = CliffordOperator(QuantumClifford.tab(Stabilizer(collect(OLD_ENTANGLE_PINS[i]))))

@testset "gate set: the old Clifford2IndexSet(:entangle)" begin
  G = DQ.old_entangle_gates()

  @testset "number of gates: exactly 19, no duplicates" begin
    @test length(G) == 19
    @test length(OLD_ENTANGLE_PINS) == 19
    @test allunique(G)
    @test allunique(pinned_tableau.(1:19))
    # Every gate is a valid Clifford (its inverse undoes it), not just 19 tableaux.
    @test all(g * inv(g) == one(CliffordOperator, 2) for g in G)
  end

  @testset "order: gate i equals pinned tableau i" begin
    # Order matters: the bond search keeps the first minimum, so it decides ties.
    @test all(G[i] == pinned_tableau(i) for i in 1:19)
    # Sensitivity: the comparison above would notice a swapped pair of gates.
    @test G[[2, 1, 3:19...]] != pinned_tableau.(1:19)
    # The pin includes the signs: flipping one changes the tableau.
    flipped = collect(OLD_ENTANGLE_PINS[17]); flipped[1] = -flipped[1]
    @test CliffordOperator(QuantumClifford.tab(Stabilizer(flipped))) != G[17]
  end

  @testset "the search matrices implement the pinned tableaux (qubit order included)" begin
    # The bond search scores `de.matrices`; the frame write uses `de.gates`. Both must be the
    # pinned gate. (Until Phase 6 the matrices were compared with the old code's ITensor
    # gates, built separately from the tableaux.) Here, independently of how they are built:
    # for M = matrices[i] acting on sites (1, 2) of a random state |a⟩ and every two-qubit
    # Pauli P, ⟨a|M† P M|a⟩ = ⟨a|C† P C|a⟩ with C = pinned tableau i (cliffconj, helpers.jl).
    # This fixes M up to a global phase, and qubit 1 of the tableau is site 1.
    de = DQ.BondSearchDisentangler(:entangle)
    s = siteinds("Qubit", 2)
    rng = Xoshiro(7)
    ops = all_pauli_ops(2)
    # ⟨P⟩ for every P in `ps` on a plain MPS (identity frame).
    expvals(mps, ps) = physical_expectations(mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(2), false), ps)
    # Largest deviation of M|a⟩ from C|a⟩ over all Paulis, for a fresh random |a⟩.
    deviation(M, C) = begin
      a = random_mps(rng, ComplexF64, s; linkdims = 2)
      b = copy(a)
      mt.apply!(b, M, [1, 2])
      maximum(abs.(expvals(b, ops) .- expvals(a, [cliffconj(inv(C), P) for P in ops])))
    end
    @test all(deviation(de.matrices[i], pinned_tableau(i)) <= 1e-12 for i in 1:19)
    @test all(de.matrices[i] ≈ Matrix(de.gates[i]) for i in 1:19)
    # Sensitivity: a different gate, or the right gate with the qubit order reversed, is
    # visibly different on a random state.
    @test all(deviation(de.matrices[i], pinned_tableau(mod1(i + 1, 19))) > 1e-3 for i in 1:19)
    rev(M) = M[[1, 3, 2, 4], [1, 3, 2, 4]]
    nonsym = [i for i in 1:19 if !(de.matrices[i] ≈ rev(de.matrices[i]))]
    @test length(nonsym) >= 10
    @test all(deviation(rev(de.matrices[i]), pinned_tableau(i)) > 1e-3 for i in nonsym)
  end

  @testset ":chi3 and :entangle search the same gates, in the same order" begin
    # In the old code `index_set` defaults to `Clifford2IndexSet(:entangle)` whatever the
    # criterion (dmrgCliffordMPS/src/Disentangler.jl), and no caller overrides it.
    for crit in (:chi3, :entangle)
      @test search_tableaux(DQ.BondSearchDisentangler(crit)) == G
      for strat in (:snake, :brickwork)
        @test search_tableaux(bond_search(DQ.sweep_disentangler(crit, strat, 3))) == G
      end
    end
  end

  @testset "argument checks" begin
    @test_throws ArgumentError DQ.BondSearchDisentangler(:full)
    @test_throws ArgumentError DQ.sweep_disentangler(:entangle, :full, 3)
    @test_throws ArgumentError DQ.sweep_disentangler(:chi2, :snake, 3)
  end
end

# -- Snake skips bonds of dimension 1 ----------------------------------------------------

# The old snake (dmrgCliffordMPS/src/Disentangler.jl) re-reads the bond dimensions at the start
# of each direction and does not search bonds of dimension 1 (a product cut has nothing to
# disentangle). A uniform-χ ≥ 2 state never has such a bond, so a snake that searches every
# bond would pass tests on those. With :entangle, the search on a dimension-1 bond finds
# tie-level "gains" (~1e-16) and applies a gate; before item 8, with no cutoff (the campssrc
# setting then), that gate also grew the bond.

"A random MPS with the bond dimensions `dims` (length N-1), behind a random frame."
function random_state_dims(rng, dims)
  N = length(dims) + 1
  s = siteinds("Qubit", N)
  mps = random_mps(rng, ComplexF64, s; linkdims = dims)
  C = random_clifford(rng, N)
  return mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(C), false)
end

@testset "snake: product MPS (all bonds of dimension 1) is left untouched" begin
  # Every bond is skipped, so no gate is applied, whatever the criterion or the cutoff:
  # 0 gates, frame bitwise unchanged, all linkdims still 1 (as in the old code).
  for crit in (:chi3, :entangle), kw in ((;), (; cutoff = 1e-12))
    @testset "($crit), $(isempty(kw) ? "no cutoff" : "cutoff 1e-12")" begin
      rng = Xoshiro(hash((:product, crit, length(kw))))
      for _ in 1:8
        N = rand(rng, 4:6)
        new = random_state_dims(rng, ones(Int, N - 1))
        frame_before = frame_dagger(new)
        num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(crit, :snake, 2; kw...), new)
        @test num == 0
        @test all(isempty(sw[3]) for sw in sweeps)
        @test gain == 0
        @test all(==(1), linkdims(mt.get_mps(new)))
        @test frame_dagger(new) == frame_before
      end
    end
  end
end

@testset "snake: mixed bond dimensions, dimension-1 bonds are skipped" begin
  # Random MPS with a run of L dimension-1 bonds (L = 1 or 3-5; the rest of dimension 2 or 3).
  # The schedule is what sees a skipped bond: a bond of dimension 1 never appears in it.
  # (Until Phase 6 this was an old-vs-new parity test; the state choice below is kept.)
  # Why runs: a skipped bond only gets a tie-level "gain" (~1e-16, a gate that keeps the cut a
  # product, e.g. SWAP) when BOTH of its sites are decoupled single-site blocks, i.e. two
  # neighbouring dimension-1 bonds on each side. A lone dimension-1 bond between entangled
  # blocks yields no such gate, so a no-skip snake is only caught by runs of L >= 3 (probe,
  # 40 states each: no-skip differs in 0/40 for L = 1, 2 and in 4-22/40 for L = 3-5).
  skipped_seen = Ref(0)
  gates = Dict{Symbol, Int}()
  for crit in (:chi3, :entangle), kw in ((;), (; cutoff = 1e-12))
    @testset "($crit), $(isempty(kw) ? "no cutoff" : "cutoff 1e-12")" begin
      rng = Xoshiro(hash((:mixed, crit, length(kw))))
      for _ in 1:8
        N = rand(rng, 8:10)
        # Run of L ones on bonds p..q. The block of sites left of the run (sites 1..p) and the
        # one to its right (q+1..N) must not have exactly 2 sites: a pure two-qubit state is
        # disentangled exactly by several gates of the set, a genuine tie decided by rounding,
        # where old-vs-new parity could not be asserted (see the file header).
        # Rejection-sample (p, q).
        local p, q
        while true
          L = rand(rng, [1, 3, 4, 5])
          p = rand(rng, 1:N-L); q = p + L - 1
          p != 2 && N - q != 2 && break
        end
        dims = rand(rng, 2:3, N - 1)
        dims[p:q] .= 1
        new = random_state_dims(rng, dims)
        d0 = linkdims(mt.get_mps(new))
        max_iter = rand(rng, 1:3)
        num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(crit, :snake, max_iter; kw...), new)
        @test 1 <= length(sweeps) <= max_iter
        @test num == sum(sw -> length(sw[3]), sweeps)
        # The rule: a gate on bond b is the only thing that changes bond b's dimension, so a
        # bond that starts at 1 is never searched and stays 1.
        searched = Set(b for sw in sweeps for (_, b) in sw[3])
        @test all(b -> d0[b] != 1, searched)
        @test all(b -> linkdims(mt.get_mps(new))[b] == 1, findall(==(1), d0))
        skipped_seen[] += count(==(1), d0)
        gates[crit] = get(gates, crit, 0) + num
      end
    end
  end
  @test skipped_seen[] >= 32     # dimension-1 bonds really were present
  @test gates[:entangle] >= 10   # not vacuous: the search acted on the other bonds
end

# -- Item 8: numerical zeros; :chi3 acts (CONVENTIONS.md §2) -----------------------------
#
# Replaces the old ":chi3 never applies a gate" testset (plan item 8: replaced, not kept
# alongside). Which mutant of the plan's list each group is meant to catch:
# - the old no-op (every value counted): every pin where :chi3 acts (Bell, cluster, two
#   Bell pairs, τ boundary, randomized-exact);
# - τ 1e6× too tight (1e-16): τ boundary at 0.1τ (gain 56 and bond 2 instead of 7 and 1);
# - τ 1e6× too loose (1e-4) and counting σ² against a σ-scale τ: τ boundary at 10τ (gain 7
#   and bond 1 instead of 56 and 2);
# - `maxdim = r` not passed (padding kept): every bond-dimension pin, `padded_bonds`;
# - a tie with the identity accepted (≤ for <): product and generic states (brickwork
#   searches every bond and would apply SWAP), and the gain pins.
# Added after the item 8 review (2026-11-20), for rules §2 rejects:
# - τ absolute (σ ≤ τ): the boundary pins on a scaled MPS;
# - τ relative to the norm √Σσ², or cumulative (search rank by discarded weight ≤ τ², or the
#   applied gate with ITensors `cutoff = τ²` instead of `maxdim = r`): the per-value testset;
# - a user `maxdim` ignored, or combined with r by max instead of min: the `maxdim` testset.

@testset "item 8 hand derivations, re-checked on the pinned tableaux (no code under test)" begin
  # Bell pair Φ⁺, stabilizers +XX, +ZZ. By hand, gate 3 (X₁ → ZY, X₂ → ZI, Z₁ → IZ, Z₂ → YZ):
  # XX → ZY·ZI = +IY and ZZ → IZ·YZ = +YI, so G₃Φ⁺ = |+i⟩|+i⟩ (Y eigenstates), a product.
  # Gate 1 (SWAP) keeps XX, ZZ; gate 2: XX → ZY·YZ = (−iX)(iX) = +XX, ZZ → IZ·ZI = +ZZ: both
  # keep it entangled. So gate 3 is the first that disentangles it; CNOT (11) is not.
  # Over all 19 (plan item 8 lists the same set, checked there on dense matrices):
  bell_set = [3, 4, 5, 7, 8, 9, 11, 13, 15, 16, 17, 18]
  @test [k for k in 1:19 if pinned_disentangles(k, P"XX", P"ZZ")] == bell_set
  @test cliffconj(pinned_tableau(3), P"XX") == P"IY" && cliffconj(pinned_tableau(3), P"ZZ") == P"YI"
  # Cluster state CZ|++⟩, stabilizers +XZ, +ZX. Gate 2: XZ → ZY·ZI = +IY, ZX → IZ·YZ = +YI,
  # a product. CNOT (11): XZ → XX·ZZ = −YY, ZX → ZI·IX = +ZX, still maximally entangled.
  cluster_set = [2, 4, 6, 7, 8, 9, 12, 13, 14, 16, 17, 18]
  @test [k for k in 1:19 if pinned_disentangles(k, P"XZ", P"ZX")] == cluster_set
  @test cliffconj(pinned_tableau(2), P"XZ") == P"IY" && cliffconj(pinned_tableau(2), P"ZX") == P"YI"
  @test cliffconj(pinned_tableau(11), P"XZ") == P"-YY"
  # Undoing a gate on product states (randomized-exact test below): gates 2–10 are undone by
  # gates 2 and 12, gates 11–19 by gates 3 and 11, SWAP by SWAP.
  @test undoers(1) == [1]
  @test all(undoers(k) == [2, 12] for k in 2:10)
  @test all(undoers(k) == [3, 11] for k in 11:19)
  # Negative checks: the helpers reject wrong answers.
  @test !pinned_disentangles(11, P"XZ", P"ZX")
  @test !local_or_swaplocal(pinned_tableau(11))                  # CNOT is not local
  @test local_or_swaplocal(tHadamard ⊗ tPhase) && local_or_swaplocal(pinned_tableau(1))
  @test !local_or_swaplocal(pinned_tableau(2) * pinned_tableau(11))   # 2 does not undo 11
  @test local_or_swaplocal(pinned_tableau(3) * pinned_tableau(11))    # 3 does
end

@testset "item 8 property helpers reject wrong answers" begin
  # A zero-padded bond: SWAP applied to a product state by `mt.apply!` (untruncated) leaves
  # bond 2 at dimension 2 holding a zero value.
  s = siteinds("Qubit", 4)
  m = MPS(ComplexF64, s, "0")
  mt.apply!(m, Matrix(mt.QubitCliffordGate(tSWAP)), [2, 3])
  @test linkdims(m) == [1, 2, 1]
  @test padded_bonds(m) == [2]
  @test isempty(padded_bonds(bell_mps(s, 2)))
  # τ is relative to σ_max, per value.
  @test numrank([1.0, 0.1TAU]) == 1 && numrank([1.0, 10TAU]) == 2
  @test numrank([2.0, 1.5TAU]) == 1 && numrank([2.0, 1.0, 0.0]) == 2
  # The spectrum helper on a known state: Φ⁺ ⊗ weak pair, bond 2 = {1, 1, ε, ε}/norm.
  σ = bond_spectrum(bell_times_weak_pair(s, 1e-3), 2)
  @test all(isapprox.(σ ./ σ[1], [1, 1, 1e-3, 1e-3]; rtol = 1e-10))
  @test abs(inner(bell_times_schmidt(s, [1, 1e-3]), bell_times_weak_pair(s, 1e-3)) - 1) <= 1e-12
  s6 = siteinds("Qubit", 6)
  σ = bond_spectrum(bell_times_schmidt(s6, [1, 0.5, 0.25, 0.125]), 3)
  @test isapprox(σ ./ σ[1], [1, 1, 0.5, 0.5, 0.25, 0.25, 0.125, 0.125]; rtol = 1e-10)
  # The rejected rules count differently from `numrank` where they should.
  @test rank_absolute([1e-6, 0.5e-10]) == 1 && rank_absolute([1e6, 2e-10]) == 2
  @test rank_vs_norm([1, 1, 1, 1, 1.5TAU]) == 4 && numrank([1, 1, 1, 1, 1.5TAU]) == 5
  @test rank_cumulative([1; fill(0.6TAU, 3)], 1) == 2           # drops two of the three
  @test rank_cumulative([1; fill(0.6TAU, 2)], 1) == 1
  # A logged gate that keeps the bond or keeps a zero fails the checks used below.
  bad = [(bond = 1, before = 2, after = 2, rank_after = 1)]
  @test !all(e.after < e.before for e in bad)
  @test !all(e.after == e.rank_after for e in bad)
end

@testset "[item 8] τ = 1e-10" begin
  @test DQ.NUMERICAL_ZERO_TAU === TAU
end

@testset "[item 8] :chi3 removes a Bell pair with the first gate of the set that does ($strat)" for
    strat in (:snake, :brickwork)
  N = 4
  s = siteinds("Qubit", N)
  ops = all_pauli_ops(N)
  ψ = mt.QubitCAMPS(bell_mps(s, 2), mt.QubitCliffordUnitary(N), false)
  @test linkdims(mt.get_mps(ψ)) == [1, 2, 1]
  ref = physical_expectations(ψ, ops)
  num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, N; min_diff = 6), ψ)
  # Bond 2: the identity has rank 2 (score 8), gate 3 rank 1 (score 1): gain 7 > min_diff,
  # so a second sweep runs, finds nothing and stops. Bonds 1 and 3 have rank 1, where no
  # gate is strictly better than the identity.
  @test num == 1 && gain == 7.0
  @test length(sweeps) == 2
  @test sweeps[1][3] == [(3, 2)] && isempty(sweeps[2][3])
  @test linkdims(mt.get_mps(ψ)) == [1, 1, 1]
  # Frame: |ψ⟩ = C|mps⟩ (§4). The MPS became G|mps⟩, so C becomes C·G† and the stored
  # C† becomes G·C† = G (C = 1 here): gate 3 on qubits (2, 3), signs included.
  @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [2, 3])
  @test frame_dagger(ψ) != embed(pinned_tableau(11), N, [2, 3])
  # The bare MPS is G₃Φ⁺ = |+i⟩|+i⟩ on sites 2, 3 (hand derivation above), |0⟩ elsewhere.
  bare = mt.QubitCAMPS(deepcopy(mt.get_mps(ψ)), mt.QubitCliffordUnitary(N), false)
  @test physical_expectations(bare, [P"IYII", P"IIYI", P"ZIII", P"IIIZ"]) ≈ [1, 1, 1, 1] atol = 1e-12
  @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= 1e-12

  # :entangle on the same state: the applied gate now also drops the zero value (the old
  # code, without cutoff, left [1, 2, 1]). Which gate wins is a rounding-decided tie.
  ψe = mt.QubitCAMPS(bell_mps(s, 2), mt.QubitCliffordUnitary(N), false)
  num_e, gain_e, _ = mt.disentangle!(DQ.sweep_disentangler(:entangle, strat, N; min_diff = 6), ψe)
  @test num_e == 1
  @test gain_e ≈ 1 atol = 1e-12
  @test linkdims(mt.get_mps(ψe)) == [1, 1, 1]
end

@testset "[item 8] :chi3 removes two Bell pairs in one sweep, gain 14 ($strat)" for
    strat in (:snake, :brickwork)
  # N = 6, Φ⁺ on qubits (1, 2) and (4, 5). Snake: bonds 1 and 4 forward (2, 3, 5 skipped).
  # Brickwork: bonds 1, 3, 5, then 4, 2. Either way gate 3 on bond 1, then on bond 4.
  N = 6
  s = siteinds("Qubit", N)
  mps = bell_mps(s, 1)
  mt.apply!(mps, Matrix(mt.QubitCliffordGate(tCNOT * (tHadamard ⊗ tId1))), [4, 5])
  ψ = mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(N), false)
  @test linkdims(mt.get_mps(ψ)) == [2, 1, 1, 2, 1]
  num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, N; min_diff = 6), ψ)
  @test num == 2 && gain == 14.0
  @test length(sweeps) == 2
  @test sweeps[1][3] == [(3, 1), (3, 4)] && isempty(sweeps[2][3])
  @test all(==(1), linkdims(mt.get_mps(ψ)))
  @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [1, 2]) * embed(pinned_tableau(3), N, [4, 5])
end

@testset "[item 8] :chi3 on a state CNOT does not disentangle: the first gate that does ($strat)" for
    strat in (:snake, :brickwork)
  # The cluster state CZ|++⟩ on qubits (2, 3) of N = 4: gate 2 is the first that makes it
  # a product (hand derivation above), CNOT is not among them.
  N = 4
  s = siteinds("Qubit", N)
  mps = MPS(ComplexF64, s, "0")
  mt.apply!(mps, Matrix(mt.QubitCliffordGate(tCPHASE * (tHadamard ⊗ tHadamard))), [2, 3])
  ψ = mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(N), false)
  @test linkdims(mt.get_mps(ψ)) == [1, 2, 1]
  # Check of the input, frame-free: the stabilizers +XZ, +ZX on qubits 2, 3.
  @test physical_expectations(ψ, [P"IXZI", P"IZXI"]) ≈ [1, 1] atol = 1e-12
  num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, N; min_diff = 6), ψ)
  @test num == 1 && gain == 7.0
  @test sweeps[1][3] == [(2, 2)] && length(sweeps) == 2 && isempty(sweeps[2][3])
  @test linkdims(mt.get_mps(ψ)) == [1, 1, 1]
  @test frame_dagger(ψ) == embed(pinned_tableau(2), N, [2, 3])
  bare = mt.QubitCAMPS(deepcopy(mt.get_mps(ψ)), mt.QubitCliffordUnitary(N), false)
  @test physical_expectations(bare, [P"IYII", P"IIYI"]) ≈ [1, 1] atol = 1e-12
end

@testset "[item 8] :chi3 on a product state: nothing applied, one empty sweep ($strat)" for
    strat in (:snake, :brickwork)
  # Every bond has rank 1, so no gate can be strictly better than the identity. (Snake
  # skips all bonds; brickwork searches them all, so a tie accepted with the identity
  # would show here.)
  rng = Xoshiro(hash((:item8_product, strat)))
  for _ in 1:4
    N = rand(rng, 4:6)
    ψ = random_state_dims(rng, ones(Int, N - 1))
    mps0 = deepcopy(mt.get_mps(ψ)); frame0 = frame_dagger(ψ)
    num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, N; min_diff = 6), ψ)
    @test num == 0 && gain == 0.0
    @test length(sweeps) == 1 && isempty(sweeps[1][3])
    @test same_tensors(mt.get_mps(ψ), mps0)
    @test frame_dagger(ψ) == frame0
  end
end

@testset "[item 8] :chi3 on a generic entangled two-qubit state: nothing applied ($strat)" for
    strat in (:snake, :brickwork)
  # A seeded random two-qubit state. No pinned gate makes it a product, checked here from
  # the tableaux alone: after gate C, the purity of qubit 1 is (1 + Σ_P ⟨C P C†⟩²)/2 over
  # P ∈ {XI, YI, ZI}, with ⟨C P C†⟩ on G|a⟩ = ⟨a|C† P C|a⟩ (`cliffconj(inv(C), P)`).
  # (cos a|00⟩ + sin a|11⟩ would not do: CNOT makes it a product.)
  rng = Xoshiro(hash((:item8_generic2, strat)))
  s = siteinds("Qubit", 2)
  for _ in 1:4
    mps = random_mps(rng, ComplexF64, s; linkdims = 2)
    bare = mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(2), false)
    purity(C) = (1 + sum(abs2, physical_expectations(bare,
      [cliffconj(inv(C), P) for P in (P"XI", P"YI", P"ZI")]))) / 2
    @test numrank(bond_spectrum(mps, 1)) == 2
    @test purity(one(CliffordOperator, 2)) < 1 - 1e-3
    @test all(purity(pinned_tableau(k)) < 1 - 1e-3 for k in 1:19)
    ψ = mt.QubitCAMPS(deepcopy(mps), mt.QubitCliffordUnitary(random_clifford(rng, 2)), false)
    mps0 = deepcopy(mt.get_mps(ψ)); frame0 = frame_dagger(ψ)
    num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, 2; min_diff = 6), ψ)
    @test num == 0 && gain == 0.0
    @test length(sweeps) == 1 && isempty(sweeps[1][3])
    @test same_tensors(mt.get_mps(ψ), mps0)
    @test frame_dagger(ψ) == frame0
  end
end

@testset "[item 8] τ boundary: a value at 0.1τ·σ_max is zero, one at 10τ·σ_max is not" begin
  # Φ⁺ on (2, 3) times (|00⟩ + ε|11⟩) on (1, 4), N = 4 (`bell_times_weak_pair`). Across
  # bond 2: identity {1, 1, ε, ε}, gates 1, 2 the same ranks (they keep Φ⁺ entangled),
  # gate 3 {1, ε} (G₃Φ⁺ is a product), all relative to σ_max. So gate 3 is chosen at both ε:
  # - ε = 0.1τ: ranks 2 → 1, scores 8 → 1, gain 7, and the 0.1τ value is dropped (bond 1);
  # - ε = 10τ:  ranks 4 → 2, scores 64 → 8, gain 56, and the 10τ value is kept (bond 2).
  # The single-bond search on bond 2 only: bonds 1 and 3 carry the weak pair too.
  N = 4
  s = siteinds("Qubit", N)
  ops = all_pauli_ops(N)
  bs = bond_search(DQ.sweep_disentangler(:chi3, :snake, N; min_diff = 6))
  for (ε, gain_exp, dim_exp, phys_tol) in ((0.1TAU, 7.0, 1, 1e-10), (10TAU, 56.0, 2, 1e-12))
    @testset "ε = $(ε / TAU)·τ" begin
      ψ = mt.QubitCAMPS(bell_times_weak_pair(s, ε), mt.QubitCliffordUnitary(N), false)
      σ0 = bond_spectrum(mt.get_mps(ψ), 2)
      @test linkdims(mt.get_mps(ψ)) == [2, 4, 2]
      @test all(isapprox.(σ0 ./ σ0[1], [1, 1, ε, ε]; rtol = 1e-6))
      ref = physical_expectations(ψ, ops)
      @test mt.disentangle!(bs, ψ, 2, false) == (1, gain_exp, 3)
      mps = mt.get_mps(ψ)
      @test linkdims(mps)[2] == dim_exp
      σ = bond_spectrum(mps, 2)
      @test numrank(σ) == dim_exp
      dim_exp == 2 && @test σ[2] / σ[1] ≈ ε rtol = 1e-6
      @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [2, 3])
      # Dropping the 0.1τ value moves expectations by O(ε); keeping it, by rounding only.
      @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= phys_tol

      # :entangle on the same state: its applied gate drops the 0.1τ value and keeps the
      # 10τ one. (Which disentangling gate wins is a rounding-decided tie.)
      ψe = mt.QubitCAMPS(bell_times_weak_pair(s, ε), mt.QubitCliffordUnitary(N), false)
      bse = bond_search(DQ.sweep_disentangler(:entangle, :snake, N; min_diff = 6))
      n_e, = mt.disentangle!(bse, ψe, 2, false)
      @test n_e == 1
      @test linkdims(mt.get_mps(ψe))[2] == dim_exp
    end
  end
end

@testset "[item 8] τ is relative to σ_max: the boundary pins on a scaled MPS (×$c)" for
    c in (1e6, 1e-6)
  # The τ-boundary state with its MPS multiplied by c: every singular value scales by c, so
  # the relative rule chooses as at c = 1 (gate 3, gains 7 and 56, bond 2 ending at 1 and 2).
  # An absolute τ (σ ≤ τ) would not: at c = 1e6 the 0.1τ values are ~7e-6 > τ and count, at
  # c = 1e-6 the 10τ values are ~7e-16 < τ and are dropped (checked on the spectra below).
  N = 4
  s = siteinds("Qubit", N)
  bs = bond_search(DQ.sweep_disentangler(:chi3, :snake, N; min_diff = 6))
  ops = all_pauli_ops(N)
  for (ε, gain_exp, dim_exp) in ((0.1TAU, 7.0, 1), (10TAU, 56.0, 2))
    @testset "ε = $(ε / TAU)·τ" begin
      mps = orthogonalize(bell_times_weak_pair(s, ε), 1)
      mps[1] *= c                                       # scale at the orthogonality centre
      ψ = mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(N), false)
      σ0 = bond_spectrum(mt.get_mps(ψ), 2)
      @test σ0[1] ≈ c / sqrt(2) rtol = 1e-10
      @test numrank(σ0) == 2 * dim_exp
      # Sensitivity: the absolute rule counts the ε values differently, at c = 1e6 for
      # ε = 0.1τ (keeps them) and at c = 1e-6 for ε = 10τ (drops them).
      (c > 1) == (ε < TAU) && @test rank_absolute(σ0) != numrank(σ0)
      ref = physical_expectations(ψ, ops) ./ c^2
      @test mt.disentangle!(bs, ψ, 2, false) == (1, gain_exp, 3)
      @test linkdims(mt.get_mps(ψ))[2] == dim_exp
      @test numrank(bond_spectrum(mt.get_mps(ψ), 2)) == dim_exp
      @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [2, 3])
      @test norm(mt.get_mps(ψ)) ≈ c rtol = 1e-10
      @test maximum(abs.(physical_expectations(ψ, ops) ./ c^2 .- ref)) <= 1e-10
    end
  end
end

@testset "[item 8] per value, relative to σ_max: not cumulative, not relative to the norm" begin
  # N = 6, `bell_times_schmidt`: Φ⁺ on (3, 4) times χ with Schmidt values μ on (1, 2 | 5, 6).
  # Bond 3 holds {1, 1} ⊗ μ; gate 3 leaves μ (the first gate that makes Φ⁺ a product, as on
  # the Bell pin). Single-bond search on bond 3. Relative to σ_max:
  # - μ = {1, a, a, a}, a = 0.6τ: identity {1, 1, a ×6}, rank 2; gate 3 {1, a, a, a}, rank 1:
  #   gain 7, bond 3 ends at 1. A cumulative rule (drop while the discarded Σσ² ≤ τ², against
  #   σ_max² or Σσ²; e.g. ITensors `cutoff = τ²` on the applied gate) keeps some of the a's.
  # - μ = {1, 1, 1, b}, b = 2τ: identity {1 ×6, b ×2}, rank 8; gate 3 {1, 1, 1, b}, rank 4:
  #   gain 512 − 64 = 448, bond 3 ends at 4. Against the norm √Σσ² ≈ 2.45·σ_max, b would be
  #   a zero on the identity (rank 6, gain 152).
  N = 6
  s = siteinds("Qubit", N)
  ops = test_paulis(Xoshiro(hash(:item8_pervalue)), N)
  a, b = 0.6TAU, 2TAU
  cases = [(name = "several values at 0.6τ", μ = [1, a, a, a], r_id = 2, r_gate = 1, gain = 7.0),
           (name = "large norm, a value at 2τ", μ = [1, 1, 1, b], r_id = 8, r_gate = 4, gain = 448.0)]
  for cs in cases
    @testset "$(cs.name)" begin
      mps = bell_times_schmidt(s, cs.μ)
      @test linkdims(mps)[3] == 8
      σ_id = bond_spectrum(mps, 3)
      @test numrank(σ_id) == cs.r_id
      # σ_gate: the spectrum gate 3 leaves (μ, normalised), from the construction.
      σ_gate = sort(cs.μ ./ norm(cs.μ); rev = true)
      @test numrank(σ_gate) == cs.r_gate
      # Sensitivity: the rejected rules give other ranks on these spectra.
      if cs.r_id == 2
        @test rank_cumulative(σ_id, σ_id[1]^2) > cs.r_id && rank_cumulative(σ_id, sum(abs2, σ_id)) > cs.r_id
        @test rank_cumulative(σ_gate, σ_gate[1]^2) > cs.r_gate && rank_cumulative(σ_gate, sum(abs2, σ_gate)) > cs.r_gate
      else
        @test rank_vs_norm(σ_id) < cs.r_id && rank_vs_norm(σ_gate) == cs.r_gate
      end
      ψ = mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(N), false)
      ref = physical_expectations(ψ, ops)
      bs = bond_search(DQ.sweep_disentangler(:chi3, :snake, N; min_diff = 6))
      @test mt.disentangle!(bs, ψ, 3, false) == (1, cs.gain, 3)
      @test linkdims(mt.get_mps(ψ))[3] == cs.r_gate
      @test numrank(bond_spectrum(mt.get_mps(ψ), 3)) == cs.r_gate
      @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [3, 4])
      # Dropping values ≤ 0.6τ·σ_max moves expectations by O(τ) at most.
      @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= 1e-9
      # :entangle: whichever gate wins (a rounding-decided tie), its SVD keeps r values.
      ψe = mt.QubitCAMPS(bell_times_schmidt(s, cs.μ), mt.QubitCliffordUnitary(N), false)
      bse = bond_search(DQ.sweep_disentangler(:entangle, :snake, N; min_diff = 6))
      @test first(mt.disentangle!(bse, ψe, 3, false)) == 1
      @test linkdims(mt.get_mps(ψe))[3] == cs.r_gate
    end
  end
end

@testset "[item 8] SVD keywords: `cutoff` truncates the applied gate only, on top of τ" begin
  # The ε = 10τ state of the τ-boundary test, with cutoff = 1e-12 (ITensors: discarded
  # weight Σσ²/Σσ²). The search ignores it: gate 3, gain 56, as without. The applied gate's
  # SVD then also drops the 10τ value, whose weight ε²/(1 + ε²) = 1e-18 is below the
  # cutoff: bond 2 ends at 1, not 2. The physical state moves by O(ε) only.
  N = 4
  s = siteinds("Qubit", N)
  ops = all_pauli_ops(N)
  ε = 10TAU
  ψ = mt.QubitCAMPS(bell_times_weak_pair(s, ε), mt.QubitCliffordUnitary(N), false)
  ref = physical_expectations(ψ, ops)
  bs = bond_search(DQ.sweep_disentangler(:chi3, :snake, N; min_diff = 6, cutoff = 1e-12))
  @test mt.disentangle!(bs, ψ, 2, false) == (1, 56.0, 3)
  @test linkdims(mt.get_mps(ψ))[2] == 1
  @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [2, 3])
  @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= 1e-8
end

@testset "[item 8] SVD keywords: a user `maxdim` caps the applied gate, combined with r by min" begin
  # §2 "SVD keywords": `maxdim` applies to the applied gate's SVD on top of the τ rule, and
  # the search ignores it. On the τ-boundary states (gate 3 chosen, as without keywords):
  # - ε = 10τ (r = 2), maxdim = 1: bond 2 ends at 1. Ignoring the user's value, or taking
  #   max(r, maxdim), would leave 2. The search still sees ranks 4 → 2: gain 56 (a search
  #   capped at 1 would see every candidate at rank 1 and apply nothing).
  # - ε = 0.1τ (r = 1), maxdim = 3: bond 2 ends at 1. Using the user's value instead of r
  #   would keep the 0.1τ value (bond 2).
  N = 4
  s = siteinds("Qubit", N)
  ops = all_pauli_ops(N)
  for (ε, maxdim, gain_exp, phys_tol) in ((10TAU, 1, 56.0, 1e-8), (0.1TAU, 3, 7.0, 1e-10))
    @testset "ε = $(ε / TAU)·τ, maxdim = $maxdim" begin
      ψ = mt.QubitCAMPS(bell_times_weak_pair(s, ε), mt.QubitCliffordUnitary(N), false)
      ref = physical_expectations(ψ, ops)
      bs = bond_search(DQ.sweep_disentangler(:chi3, :snake, N; min_diff = 6, maxdim))
      @test mt.disentangle!(bs, ψ, 2, false) == (1, gain_exp, 3)
      @test linkdims(mt.get_mps(ψ))[2] == 1
      @test frame_dagger(ψ) == embed(pinned_tableau(3), N, [2, 3])
      @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= phys_tol
    end
  end
end

@testset "[item 8] randomized-exact: Clifford-structured states go back to a product ($strat)" for
    strat in (:snake, :brickwork)
  # A random product state, a random gate k ∈ 2:19 of the set on each of some disjoint qubit
  # pairs, a random frame. Each pair has rank 2, and the first gate that makes it a product
  # again is first(undoers(k)) (2 or 3; derived above from the tableaux). One sweep does all
  # pairs, gaining 7 each; the second finds nothing. The physical state is unchanged (all 4^N
  # Pauli strings up to N = 5, 300 random ones beyond).
  rng = Xoshiro(hash((:item8_structured, strat)))
  for _ in 1:6
    N = rand(rng, 4:8)
    bonds = disjoint_bonds(rng, N)
    ks = rand(rng, 2:19, length(bonds))
    ψ = structured_state(rng, N, bonds, ks)
    @test linkdims(mt.get_mps(ψ)) == [b in bonds ? 2 : 1 for b in 1:N-1]
    ops = test_paulis(rng, N)
    ref = physical_expectations(ψ, ops)
    num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, N; min_diff = 6), ψ)
    @test all(==(1), linkdims(mt.get_mps(ψ)))
    @test num == length(bonds) && gain == 7.0 * length(bonds)
    @test length(sweeps) == 2 && isempty(sweeps[2][3])
    # Order of the sweep: snake goes forward; brickwork does odd bonds, then even in reverse.
    order = strat == :snake ? sort(bonds) : [sort(filter(isodd, bonds)); sort(filter(iseven, bonds); rev = true)]
    expected = [(first(undoers(ks[findfirst(==(b), bonds)])), b) for b in order]
    @test sweeps[1][3] == expected
    @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= 1e-12
  end
end

@testset "[item 8] properties at every applied gate and after every sweep ($crit, $strat)" for
    crit in (:chi3, :entangle), strat in (:snake, :brickwork)
  # On inputs without numerical zeros, at every gate the bond search applies:
  # - the bond keeps exactly its values above τ·σ_max (stored dimension = numerical rank);
  # - for :chi3, its dimension strictly drops (r_gate < r_identity ≤ dimension before).
  # After every sweep: no bond holds a numerical zero, and the physical state is unchanged.
  # States: Clifford-structured (as above) and generic (random bond dimensions 2–3), N = 4–6,
  # one sweep at a time, three sweeps each.
  rng = Xoshiro(hash((:item8_properties, crit, strat)))
  gates_seen = 0
  for trial in 1:6
    N = rand(rng, 4:6)
    ψ = if isodd(trial)
      bonds = disjoint_bonds(rng, N)
      structured_state(rng, N, bonds, rand(rng, 2:19, length(bonds)))
    else
      # Capped at the largest possible rank (2 on the end bonds), so there is no padding.
      random_state_dims(rng, [min(rand(rng, 2:3), 2^b, 2^(N - b)) for b in 1:N-1])
    end
    @test isempty(padded_bonds(mt.get_mps(ψ)))       # the input has no numerical zeros
    ops = test_paulis(rng, N)
    ref = physical_expectations(ψ, ops)
    for _ in 1:3
      de, log = recording_sweep(crit, strat, 1; min_diff = -1.0)
      # The logging wrapper is transparent: same result as the plain sweep on a copy (same
      # gates, frame and state; new link indices get fresh ids, so compare by overlap).
      ψc = deepcopy(ψ)
      num_c, gain_c, info_c = mt.disentangle!(DQ.sweep_disentangler(crit, strat, 1; min_diff = -1.0), ψc)
      num, gain, info = mt.disentangle!(de, ψ)
      @test (num, gain, info) == (num_c, gain_c, info_c)
      @test frame_dagger(ψ) == frame_dagger(ψc)
      @test linkdims(mt.get_mps(ψ)) == linkdims(mt.get_mps(ψc))
      @test abs(inner(mt.get_mps(ψ), mt.get_mps(ψc)) - 1) <= 1e-12
      @test length(log) == num
      gates_seen += num
      @test all(e.after == e.rank_after for e in log)
      crit == :chi3 && @test all(e.after < e.before for e in log)
      @test isempty(padded_bonds(mt.get_mps(ψ)))
      @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= 1e-12
    end
  end
  @test gates_seen >= 3                               # not vacuous
end

# -- Physical invariance -----------------------------------------------------------------

@testset "every sweep leaves the physical state unchanged" begin
  # The disentangler rewrites |ψ⟩ = C|mps⟩ as (C·G†)(G|mps⟩): every Pauli expectation value
  # of the physical state must be unchanged after every single sweep (only numerical zeros,
  # ≤ τ·σ_max, are discarded: rounding level). N = 4: all 256 strings. Generic states and Clifford-structured states
  # (where :entangle's choices are tie-sensitive, but the physics is not).
  N = 4
  ops = all_pauli_ops(N)
  s = siteinds("Qubit", N)
  rng = Xoshiro(2026)
  structured() = begin
    # Bell pairs on bonds 1 and 3 behind a random frame: Clifford-structured, with ties.
    ψ = bell_mps(s, 1)
    mt.apply!(ψ, Matrix(mt.QubitCliffordGate(tCNOT * (tHadamard ⊗ tId1))), [3, 4])
    mt.QubitCAMPS(ψ, mt.QubitCliffordUnitary(random_clifford(rng, N)), false)
  end
  states = [(name = "generic", make = () -> random_state(rng, N, 3)),
            (name = "Clifford-structured", make = structured)]
  for st in states, crit in (:entangle, :chi3), strat in (:snake, :brickwork)
    @testset "$(st.name), ($crit, $strat)" begin
      acted = 0
      for _ in 1:2
        ψ = st.make()
        ref = physical_expectations(ψ, ops)
        @test maximum(abs, ref) > 0.1                    # not a trivial state
        de = DQ.sweep_disentangler(crit, strat, 1; min_diff = -1.0)  # exactly one sweep
        for _ in 1:3
          num, _, _ = mt.disentangle!(de, ψ)
          acted += num
          @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) <= 1e-12
        end
        @test abs(norm(mt.get_mps(ψ)) - 1) <= 1e-12
      end
      crit == :entangle && @test acted >= 1              # the frame really changed
      # [item 8] :chi3 now acts on the two Bell pairs (each rank 2 → 1). The old
      # `crit == :chi3 && @test acted == 0` (the no-op) is replaced, not kept.
      if crit == :chi3 && st.name == "Clifford-structured"
        @test acted >= 1
      end
    end
  end
  # Sensitivity: applying a gate to the MPS without updating the frame is detected.
  ψ = random_state(rng, N, 3)
  ref = physical_expectations(ψ, ops)
  orthogonalize!(mt.get_mps(ψ), 2)
  mt.apply!(mt.get_mps(ψ), Matrix(mt.QubitCliffordGate(tCNOT)), [2, 3])
  @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) > 1e-3
end
