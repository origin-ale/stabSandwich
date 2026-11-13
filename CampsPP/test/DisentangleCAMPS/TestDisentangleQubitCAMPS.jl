# Layer-end disentangler for MagicTensors `QubitCAMPS` (CampsPP/src/DisentangleQubitCAMPS.jl,
# plan Phase 3). The specification here is the old CliffordMPS `disentangle!`: the port must
# make the same choices (plan, "Phase 3 status"). Expected values come from the old code while
# it exists, from hand-checked tableau pins, and from frame-invariant physics; never from the
# new implementation.
#
# Decisions this file encodes (user, 2026-11-13):
# - `:chi3` stays the faithful no-op of the old code (plan, "Note: the old `:chi3`
#   disentangler never acts"). It is asserted here as a fact about the algorithm. If `:chi3` is
#   ever made meaningful, the no-op testset below must be replaced, not loosened.
# - `:entangle` keeps the old tie rule (no ε). On Clifford-structured states the chosen gate
#   depends on rounding, so old-vs-new parity of frames and schedules is asserted only on
#   generic random states; Clifford-structured states are checked through frame-invariant
#   quantities.

using ITensors, ITensorMPS
using Random
using LinearAlgebra: diag, norm
import CampsPP
import CliffordMPS as cmps
import MagicTensors as mt

const DQ = CampsPP.DisentangleCAMPS

# -- Helpers -----------------------------------------------------------------------------

"Gate i (1-based) of the old `Clifford2IndexSet(:entangle)`, as an old index and a tableau."
old_entangle_index(i) = cmps.Clifford2Index(cmps.Clifford2IndexSet(:entangle), i - 1)
old_entangle_tableau(i) = cmps.clifford_gate(old_entangle_index(i))

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

"The stored frame C† of a CAMPS, as a copy (old: the `Cdag` field; new: fetched daggered)."
frame_dagger(ψ::cmps.CAMPS) = copy(ψ.Cdag)
frame_dagger(ψ::mt.QubitCAMPS) = copy(mt.get_clifford(ψ; daggered = true).clifford)

"""
The same generic state in both stacks: a random MPS (independent copies sharing site
indices) behind a random frame C.
"""
function random_state_pair(rng, N, χ)
  s = siteinds("Qubit", N)
  mps = random_mps(rng, ComplexF64, s; linkdims = χ)
  C = random_clifford(rng, N)
  old = cmps.CAMPS(copy(mps), inv(C))
  new = mt.QubitCAMPS(copy(mps), mt.QubitCliffordUnitary(copy(C)), false)
  return old, new
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
    @test length(cmps.Clifford2IndexSet(:entangle)) == 19
    @test length(G) == 19
    @test length(OLD_ENTANGLE_PINS) == 19
    @test allunique(G)
    @test allunique(pinned_tableau.(1:19))
    # Every gate is a valid Clifford (its inverse undoes it), not just 19 tableaux.
    @test all(g * inv(g) == one(CliffordOperator, 2) for g in G)
  end

  @testset "order: gate i equals old gate i (old code)" begin
    # Order matters: the bond search keeps the first minimum, so it decides ties.
    @test all(G[i] == old_entangle_tableau(i) for i in 1:19)
    # The old index numbers, in order (schedules of the old code carry these numbers).
    @test [cmps.clifford_number(old_entangle_index(i)) for i in 1:19] == cmps.entangle_set
    # Sensitivity: the comparison above would notice a swapped pair of gates.
    @test G[[2, 1, 3:19...]] != [old_entangle_tableau(i) for i in 1:19]
  end

  @testset "order: gate i equals pinned tableau i" begin
    @test all(G[i] == pinned_tableau(i) for i in 1:19)
    @test all(old_entangle_tableau(i) == pinned_tableau(i) for i in 1:19)
    # The pin includes the signs: flipping one changes the tableau.
    flipped = collect(OLD_ENTANGLE_PINS[17]); flipped[1] = -flipped[1]
    @test CliffordOperator(QuantumClifford.tab(Stabilizer(flipped))) != G[17]
  end

  @testset "the search matrices are the old MPS gates (up to a global phase)" begin
    # The bond search scores `de.matrices`; the frame write uses `de.gates`. Both must be the
    # old gate, and in the old code the MPS side used ITensor gates built separately from the
    # tableaux, so compare their action on a random two-qubit state (qubit order included).
    de = DQ.BondSearchDisentangler(:entangle)
    s = siteinds("Qubit", 2)
    rng = Xoshiro(7)
    overlap(i, j) = begin
      ψ = random_mps(rng, ComplexF64, s; linkdims = 2)
      a = ITensors.apply(cmps.clifford_gate(old_entangle_index(i), s[1], s[2]), copy(ψ))
      b = copy(ψ)
      mt.apply!(b, de.matrices[j], [1, 2])
      abs(inner(a, b))
    end
    @test all(abs(1 - overlap(i, i)) <= 1e-12 for i in 1:19)
    @test all(de.matrices[i] ≈ Matrix(de.gates[i]) for i in 1:19)
    # Sensitivity: a different gate, or the right gate with the qubit order reversed, is
    # visibly different on a random state.
    @test all(overlap(i, mod1(i + 1, 19)) < 1 - 1e-3 for i in 1:19)
    rev(M) = M[[1, 3, 2, 4], [1, 3, 2, 4]]
    nonsym = [i for i in 1:19 if !(de.matrices[i] ≈ rev(de.matrices[i]))]
    @test length(nonsym) >= 10
    @test all(begin
      ψ = random_mps(rng, ComplexF64, s; linkdims = 2)
      a = ITensors.apply(cmps.clifford_gate(old_entangle_index(i), s[1], s[2]), copy(ψ))
      b = copy(ψ); mt.apply!(b, rev(de.matrices[i]), [1, 2])
      abs(inner(a, b)) < 1 - 1e-3
    end for i in nonsym)
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

# -- Old-vs-new parity on generic states -------------------------------------------------

@testset "old vs new disentangle! on generic states" begin
  # Generic random MPS behind random frames: no rounding ties, so both stacks must make the
  # same choices (Phase 3 probe: 200/200 identical). Compared: number of sweeps, gain per
  # sweep, schedule (old index number and bond of every applied gate, in order), the frame
  # (bitwise), the bond dimensions and the MPS.
  gates_applied = Dict{Symbol, Int}()
  for crit in (:chi3, :entangle), strat in (:snake, :brickwork),
      kw in ((;), (; cutoff = 1e-12)), min_diff in (0.0, 1e-6, 6.0)
    rng = Xoshiro(hash((crit, strat, length(kw), min_diff)))
    @testset "($crit, $strat), $(isempty(kw) ? "no cutoff" : "cutoff 1e-12"), min_diff = $min_diff" begin
      for _ in 1:3
        N = rand(rng, 4:7); χ = rand(rng, 2:4); max_iter = rand(rng, 1:4)
        old, new = random_state_pair(rng, N, χ)
        diffs, schedules, _ = cmps.disentangle!(old, cmps.DisentangleStrategy(strat), max_iter;
          criterion = cmps.DisentangleCriterion(crit), min_diff, kw...)
        num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(crit, strat, max_iter; min_diff, kw...), new)

        @test length(sweeps) == length(diffs) <= max_iter
        @test all(abs(d - sw[2]) <= 1e-12 for (d, sw) in zip(diffs, sweeps))
        new_schedules = [[(cmps.entangle_set[g], b, b + 1) for (g, b) in sw[3]] for sw in sweeps]
        @test [Tuple.(sch) for sch in schedules] == new_schedules
        @test num == sum(length, schedules)
        @test gain ≈ sum(diffs) atol = 1e-12
        @test frame_dagger(old) == frame_dagger(new)
        @test linkdims(old.mps) == linkdims(mt.get_mps(new))
        @test abs(inner(old.mps, mt.get_mps(new))) ≈ 1 atol = 1e-12
        gates_applied[crit] = get(gates_applied, crit, 0) + num
      end
    end
  end
  # Not vacuous: :entangle did act on these states (and :chi3 did not, see below).
  @test gates_applied[:entangle] >= 50
  @test gates_applied[:chi3] == 0
end

# -- :chi3 is a no-op --------------------------------------------------------------------

@testset ":chi3 never applies a gate (faithful to the old code)" begin
  # Mechanism (plan, "Note: the old `:chi3` disentangler never acts"): the search SVDs are
  # untruncated, and ITensors keeps the zero singular values, so every candidate gate scores
  # (min(2χ_L, 2χ_R))³, the same as the identity, and none is strictly better.
  N = 4
  s = siteinds("Qubit", N)

  @testset "mechanism: the untruncated search SVD keeps zero singular values" begin
    # CNOT turns the Bell pair on bond 2 into a product state: one nonzero singular value,
    # but the search still counts two.
    m = bell_mps(s, 2)
    orthogonalize!(m, 2)
    _, S, _ = mt.apply_and_svd(m, Matrix(mt.QubitCliffordGate(tCNOT)), 2)
    σ = diag(S)
    @test length(σ) == 2
    @test maximum(σ) ≈ 1 atol = 1e-12
    @test minimum(σ) <= 1e-14
  end

  @testset "a Bell pair that one gate of the set removes: ($strat), $(isempty(kw) ? "no cutoff" : "cutoff 1e-12")" for
      strat in (:snake, :brickwork), kw in ((;), (; cutoff = 1e-12))
    # :chi3 leaves the state untouched, although CNOT (gate 11) would cut bond 2 from 2 to 1.
    ψ = mt.QubitCAMPS(bell_mps(s, 2), mt.QubitCliffordUnitary(N), false)
    mps0 = deepcopy(mt.get_mps(ψ)); frame0 = frame_dagger(ψ)
    num, gain, sweeps = mt.disentangle!(DQ.sweep_disentangler(:chi3, strat, N; kw...), ψ)
    @test num == 0 && gain == 0.0
    @test length(sweeps) == 1 && isempty(sweeps[1][3])     # stops after one empty sweep
    @test same_tensors(mt.get_mps(ψ), mps0)
    @test frame_dagger(ψ) == frame0

    # The old code does the same.
    old = cmps.CAMPS(bell_mps(s, 2))
    diffs, schedules, _ = cmps.disentangle!(old, cmps.DisentangleStrategy(strat), N;
      criterion = cmps.DisentangleCriterion(:chi3), kw...)
    @test diffs == [0.0] && all(isempty, schedules)
    @test old.Cdag == one(CliffordOperator, N)

    # Control: :entangle on the same state does remove the pair (gain 1 bit). With the
    # cutoff, the applied gate's SVD also drops the zero singular value.
    ψe = mt.QubitCAMPS(bell_mps(s, 2), mt.QubitCliffordUnitary(N), false)
    num_e, gain_e, _ = mt.disentangle!(DQ.sweep_disentangler(:entangle, strat, N; kw...), ψe)
    @test num_e >= 1
    @test gain_e ≈ 1 atol = 1e-12
    @test linkdims(mt.get_mps(ψe)) == (isempty(kw) ? [1, 2, 1] : [1, 1, 1])
  end

  @testset "random generic states, any max_iter and min_diff: ($strat)" for strat in (:snake, :brickwork)
    rng = Xoshiro(hash((:chi3_noop, strat)))
    for _ in 1:4
      N = rand(rng, 3:6)
      _, ψ = random_state_pair(rng, N, rand(rng, 2:4))
      mps0 = deepcopy(mt.get_mps(ψ)); frame0 = frame_dagger(ψ)
      num, gain, sweeps = mt.disentangle!(
        DQ.sweep_disentangler(:chi3, strat, 5; min_diff = -1.0, cutoff = 1e-12), ψ)
      @test num == 0 && gain == 0.0
      @test length(sweeps) == 5                            # min_diff < 0: every sweep runs
      @test same_tensors(mt.get_mps(ψ), mps0)
      @test frame_dagger(ψ) == frame0
    end
  end
end

# -- Physical invariance -----------------------------------------------------------------

@testset "every sweep leaves the physical state unchanged" begin
  # The disentangler rewrites |ψ⟩ = C|mps⟩ as (C·G†)(G|mps⟩): every Pauli expectation value
  # of the physical state must be unchanged after every single sweep (no cutoff, so nothing
  # is discarded). N = 4: all 256 strings. Generic states and Clifford-structured states
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
  states = [(name = "generic", make = () -> random_state_pair(rng, N, 3)[2]),
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
      crit == :chi3 && @test acted == 0
    end
  end
  # Sensitivity: applying a gate to the MPS without updating the frame is detected.
  ψ = random_state_pair(rng, N, 3)[2]
  ref = physical_expectations(ψ, ops)
  orthogonalize!(mt.get_mps(ψ), 2)
  mt.apply!(mt.get_mps(ψ), Matrix(mt.QubitCliffordGate(tCNOT)), [2, 3])
  @test maximum(abs.(physical_expectations(ψ, ops) .- ref)) > 1e-3
end
