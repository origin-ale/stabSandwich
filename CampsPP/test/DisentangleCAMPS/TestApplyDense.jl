# Randomized-exact tests of `apply!` and `evolve*` against the dense reference (T4/T6).
# Seeded random circuits; each instance is compared exactly, through all 4^N Pauli
# expectation values, so the comparison is frame-invariant (CONVENTIONS.md §4).
#
# Expected values come from DenseReference.jl only. The code's own `paulinature` is used
# for one thing: counting which branches the random circuits reached.
#
# Every test runs on both stacks (STATE_TYPES, helpers.jl): the old CliffordMPS `CAMPS` and
# the MagicTensors `QubitCAMPS`. The hard-spot-6 special-angle cases are `@test_broken` on the
# old type only (`broken_on_old`) and must pass on `QubitCAMPS`.

using CampsPP.DisentangleCAMPS
using CampsPP.DisentangleCAMPS: paulinature, evolve_bonddim
using CliffordMPS
using QuantumClifford
using LinearAlgebra: norm
using Random
using Test

isdefined(@__MODULE__, :DenseReference) || include("DenseReference.jl")
const DR = DenseReference

const EXACT_ATOL = 1e-10

# -- Conversions ------------------------------------------------------------------------

"Unsigned string over IXYZ → Hermitian PauliOperator (phase byte 0)."
pauliop(str::AbstractString) =
  PauliOperator(0x0, [c in ('X', 'Y') for c in str], [c in ('Z', 'Y') for c in str])

"PauliOperator with phase ±1 → signed string for the dense reference."
function paulistr(P::PauliOperator)
  P.phase[] in (0x0, 0x2) || throw(ArgumentError("non-Hermitian Pauli"))
  letters = [x ? (z ? 'Y' : 'X') : (z ? 'Z' : 'I') for (x, z) in zip(xbit(P), zbit(P))]
  return (P.phase[] == 0x2 ? "-" : "") * String(letters)
end

@testset "assumption pins: Pauli ↔ string" begin
  @test pauliop("IXYZ") == P"IXYZ"
  @test pauliop("Y") == P"Y"
  @test paulistr(P"-XIZY") == "-XIZY"
  @test paulistr(P"YZ") == "YZ"
end

# -- Random circuits --------------------------------------------------------------------

const CLIFFORD_ANGLES = [k * π / 4 for k in -4:4]

"A random angle in (−π, π), at least 0.05 away from every Clifford angle."
function general_angle(rng)
  while true
    ϕ = (2rand(rng) - 1) * π
    all(abs(ϕ - c) > 0.05 for c in CLIFFORD_ANGLES) && return ϕ
  end
end

"A random non-identity N-qubit Pauli string with a random sign."
function random_gate_pauli(rng, N)
  while true
    letters = rand(rng, ['I', 'X', 'Y', 'Z'], N)
    all(==('I'), letters) && continue
    P = pauliop(String(letters))
    return rand(rng, Bool) ? P : -P
  end
end

"`depth` random gates; each angle is drawn from `angle_sources` (numbers or `:general`)."
function random_circuit(rng, N, depth, angle_sources)
  paulis = [random_gate_pauli(rng, N) for _ in 1:depth]
  phis = map(1:depth) do _
    src = rand(rng, angle_sources)
    src === :general ? general_angle(rng) : Float64(src)
  end
  return paulis, phis
end

# -- Comparison -------------------------------------------------------------------------

"Largest |⟨P⟩_CAMPS − ⟨P⟩_dense| over all 4^N strings; also checks ⟨P⟩_CAMPS is real."
function max_deviation(ψ, dense)
  N = nqubits_state(ψ)
  dev = 0.0
  for (str, ref) in zip(DR.all_paulis(N), DR.all_expectations(dense))
    ev = state_expectation(ψ, pauliop(str))
    dev = max(dev, abs(real(ev) - ref), abs(imag(ev)))
  end
  return dev
end

dense_state(N, paulis, phis) = DR.run_circuit(DR.zero_state(N), paulistr.(paulis), phis)

"""
Run the 4-argument `apply!` over the circuit from |0…0⟩, on a state of type T.
Returns (ψ, k, natures, bookkeeping_ok, n_trivial_with_k>0).

Frame-dependent bookkeeping, checked gate by gate against CONVENTIONS.md §4 with the
independent `spec_nature`: Q = C†PC (`heisenberg_image`, helpers.jl), `paulinature` must
agree with the spec rule, and k must grow by one exactly for a `:disentanglable` gate at an
angle that T does not take as Clifford (`recognised_clifford`, helpers.jl).
"""
function run_camps(T, N, paulis, phis)
  ψ = fresh_state(T, N)
  k = 0
  kspec = 0
  natures = Symbol[]
  ok = true
  ntriv = 0
  for (P, ϕ) in zip(paulis, phis)
    nature = spec_nature(k, heisenberg_image(ψ, P))
    push!(natures, nature)
    ntriv += nature == :trivial && k > 0
    ok &= paulinature(k, inv(frame_dagger(ψ)), P) == nature
    kspec += nature == :disentanglable && !recognised_clifford(T, ϕ)
    k = apply!(ψ, k, P, ϕ)
    ok &= k == kspec
  end
  return ψ, k, natures, ok, ntriv
end

"The truncating `apply!` with a negligible cutoff over the circuit; same return shape as `run_camps`."
function run_truncating(T, N, paulis, phis)
  ψ = fresh_state(T, N)
  for (P, ϕ) in zip(paulis, phis)
    apply!(ψ, P, ϕ, 1e-15)
  end
  return ψ, 0, Symbol[], true, 0
end

"Seed for a named test set: stable across Julia versions, different for every name."
setseed(base, name) = Xoshiro(base + sum(Int, codeunits(name)))

"Worst deviation over `ncirc` random circuits, plus the fraction of circuits within tolerance."
function sweep(T, rng, N, depth, ncirc, angle_sources; runner = run_camps)
  worst = 0.0
  nok = 0
  natures = Symbol[]
  bookkeeping = true
  ntriv = 0
  for _ in 1:ncirc
    paulis, phis = random_circuit(rng, N, depth, angle_sources)
    ψ, k, nat, ok, nt = runner(T, N, paulis, phis)
    append!(natures, nat)
    bookkeeping &= ok && 0 <= k <= N
    ntriv += nt
    dev = max_deviation(ψ, dense_state(N, paulis, phis))
    worst = max(worst, dev)
    nok += dev <= EXACT_ATOL
  end
  return (worst = worst, nok = nok, n = ncirc, natures = natures,
          bookkeeping = bookkeeping, ntrivial_kpos = ntriv)
end

# -- Tests ------------------------------------------------------------------------------

@testset "apply!(ψ, k, P, ϕ) vs dense, $(stackname(T))" for T in STATE_TYPES
  @testset "general angles" begin
    rng = Xoshiro(1101)
    natures = Symbol[]
    ntrivial_kpos = 0
    for (N, depth, ncirc) in [(2, 8, 40), (3, 10, 60), (4, 12, 10)]
      r = sweep(T, rng, N, depth, ncirc, [:general])
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
      @test r.bookkeeping                        # frame-dependent: natures and k
      append!(natures, r.natures)
      ntrivial_kpos += r.ntrivial_kpos
    end
    # Every `paulinature` branch was reached (no vacuous coverage).
    @test count(==(:disentanglable), natures) >= 50
    @test count(==(:logical), natures) >= 50
    @test count(==(:trivial), natures) >= 10
    @test ntrivial_kpos >= 5                     # `:trivial` with magic qubits present
  end

  @testset "k reaches N" begin
    # After enough general rotations every qubit is a magic qubit.
    rng = Xoshiro(1102)
    paulis, phis = random_circuit(rng, 3, 12, [:general])
    ψ, k, _, ok, _ = run_camps(T, 3, paulis, phis)
    @test k == 3
    @test ok
    @test max_deviation(ψ, dense_state(3, paulis, phis)) <= EXACT_ATOL
  end

  # Clifford branches that agree with exp(iϕP) in both stacks, alone and mixed with
  # general angles.
  @testset "correct Clifford branches: $name" for (name, sources) in [
    ("+π/4 only", [π/4]),
    ("−3π/4 only", [-3π/4]),
    ("0 and ±π only", [0.0, π, -π]),
    ("±π/2 only", [π/2, -π/2]),
    ("general with +π/4", [:general, π/4]),
    ("general with −3π/4", [:general, -3π/4]),
    ("general with 0, ±π", [:general, 0.0, π, -π]),
    ("general with +π/4, −3π/4, 0, π", [:general, :general, π/4, -3π/4, 0.0, π]),
  ]
    rng = setseed(1103, name)
    for (N, depth, ncirc) in [(2, 8, 20), (3, 10, 20)]
      r = sweep(T, rng, N, depth, ncirc, sources)
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
      @test r.bookkeeping                        # frame-dependent: natures and k
    end
  end

  # Known bugs in the old `apply!` (plan, hard spot 6; CONVENTIONS.md §2): `@test_broken` on
  # the old type, fixed in `QubitCAMPS` (Phase 2).
  @testset "special-angle cases (hard spot 6): $name" for (name, sources) in [
    ("−π/4 only", [-π/4]),
    ("+3π/4 only", [3π/4]),
    ("±π/4 mixed", [π/4, -π/4]),
    ("general with −π/4", [:general, -π/4]),
    ("general with +3π/4", [:general, 3π/4]),
    ("general with ±π/2", [:general, π/2, -π/2]),
    ("+π/4 with ±π/2", [π/4, π/2, -π/2]),
  ]
    rng = setseed(1104, name)
    r = sweep(T, rng, 3, 10, 20, sources)
    if broken_on_old(T)
      @test_broken r.nok == r.n
    else
      @test r.nok == r.n
      @test r.worst <= EXACT_ATOL
      @test r.bookkeeping                        # frame-dependent: natures and k
    end
  end
end

@testset "apply!(ψ, P, ϕ, thl) vs dense (no truncation), $(stackname(T))" for T in STATE_TYPES
  # The truncating method with a negligible cutoff must be exact too.
  @testset "$name" for (name, sources) in [
    ("general", [:general]),
    ("general with +π/4, −3π/4, 0, π", [:general, :general, π/4, -3π/4, 0.0, π]),
    # Without π/4 gates this method never changes the frame except by Paulis, and the old
    # ±π/2 branch is then correct (Paulis commute as Clifford operators).
    ("general with ±π/2", [:general, π/2, -π/2]),
  ]
    rng = setseed(1105, name)
    for (N, depth, ncirc) in [(2, 8, 20), (3, 10, 20), (4, 12, 5)]
      r = sweep(T, rng, N, depth, ncirc, sources; runner = run_truncating)
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
    end
  end

  @testset "special-angle cases (hard spot 6): $name" for (name, sources) in [
    ("general with −π/4", [:general, -π/4]),
    ("general with +3π/4", [:general, 3π/4]),
    ("+π/4 with ±π/2", [π/4, π/2, -π/2]),
    ("general with +π/4 and ±π/2", [:general, π/4, π/2, -π/2]),
  ]
    rng = setseed(1106, name)
    r = sweep(T, rng, 3, 10, 20, sources; runner = run_truncating)
    if broken_on_old(T)
      @test_broken r.nok == r.n
    else
      @test r.nok == r.n
      @test r.worst <= EXACT_ATOL
    end
  end
end

@testset "evolve and evolve_bonddim vs dense, $(stackname(T))" for T in STATE_TYPES
  bonddim = DisentangleCAMPS.bonddim
  rng = Xoshiro(1107)
  N = 3
  paulis, phis = random_circuit(rng, N, 12, [:general, :general, π/4])

  # evolve(ψ, t, …) applies the first t gates and returns the same k as gate-by-gate `apply!`.
  for t in (0, 5, 12)
    ψ, k = evolve(fresh_state(T, N), t, paulis, phis)
    @test max_deviation(ψ, dense_state(N, paulis[1:t], phis[1:t])) <= EXACT_ATOL
    _, kref, _, ok, _ = run_camps(T, N, paulis[1:t], phis[1:t])
    @test ok && k == kref
  end

  # evolve_bonddim with χ above the maximum possible bond dimension (2 for N = 3) runs
  # the whole circuit and reports s = number of gates applied.
  ψ, k, s = evolve_bonddim(fresh_state(T, N), 64, paulis, phis)
  @test s == length(paulis)
  @test k == run_camps(T, N, paulis, phis)[2]
  @test max_deviation(ψ, dense_state(N, paulis, phis)) <= EXACT_ATOL

  # With χ = 2 it stops at the first gate after which the bond dimension is ≥ χ
  # (frame-dependent): reached after s gates, not yet after s − 1.
  ψ, k, s = evolve_bonddim(fresh_state(T, N), 2, paulis, phis)
  @test 1 < s < length(paulis)
  @test bonddim(ψ) >= 2
  @test bonddim(evolve(fresh_state(T, N), s - 1, paulis, phis)[1]) < 2
  @test k == run_camps(T, N, paulis[1:s], phis[1:s])[2]
  @test max_deviation(ψ, dense_state(N, paulis[1:s], phis[1:s])) <= EXACT_ATOL
end

@testset "evolve_bonddim: hand-built stopping point, $(stackname(T))" for T in STATE_TYPES
  # N = 2 from |00⟩. Gates 1 and 2 are single-qubit X rotations: each is `:disentanglable`
  # with trivial disentangler (no other qubit involved), so the frame stays the identity,
  # the MPS is the product exp(i·0.4·X)|0⟩ ⊗ exp(i·0.7·X)|0⟩ (bond dimension 1) and k = 2.
  # Gate 3, exp(i·0.5·ZZ), is `:logical` and entangles the two qubits: bond dimension 2.
  # (ZZ is not a product operator on that state: ⟨ZZ⟩ = cos 0.8 · cos 1.4 ≠ ±1.)
  bonddim = DisentangleCAMPS.bonddim
  paulis = [P"XI", P"IX", P"ZZ", P"XY"]
  phis = [0.4, 0.7, 0.5, 0.3]

  ψ, k, s = evolve_bonddim(fresh_state(T, 2), 2, paulis, phis)
  @test s == 3                                   # stops right after the entangling gate
  @test k == 2
  @test bonddim(ψ) == 2
  @test max_deviation(ψ, dense_state(2, paulis[1:3], phis[1:3])) <= EXACT_ATOL

  ψ, k, s = evolve_bonddim(fresh_state(T, 2), 3, paulis, phis)   # χ = 3 cannot be reached at N = 2
  @test s == 4
  @test k == 2
  @test max_deviation(ψ, dense_state(2, paulis, phis)) <= EXACT_ATOL

  ψ, k, s = evolve_bonddim(fresh_state(T, 2), 1, paulis, phis)   # already at χ: nothing is applied
  @test s == 0 && k == 0
end

@testset "angles next to Clifford angles, $(stackname(T))" for T in STATE_TYPES
  # CONVENTIONS.md §2: the state must be exp(iϕP)|ψ⟩ for every real ϕ. An implementation
  # that rounds ϕ = c ± 1e-3 to the Clifford angle c is wrong by about 2e-3.
  multiples = [m * π / 4 for m in -4:4]

  @testset "N = 1: state and k" begin
    for c in multiples, δ in (-1e-3, 1e-3)
      ϕ = c + δ
      ψ = fresh_state(T, 1)
      k = apply!(ψ, 0, P"X", ϕ)
      @test k == 1                               # frame-dependent: not treated as Clifford
      @test max_deviation(ψ, dense_state(1, [P"X"], [ϕ])) <= EXACT_ATOL
    end
  end

  @testset "N = 2: inside a circuit" begin
    worst = 0.0
    for c in multiples, δ in (-1e-3, 1e-3)
      paulis = [P"XY", P"ZX", P"-YI"]
      phis = [0.37, c + δ, 0.9]
      ψ, _, _, ok, _ = run_camps(T, 2, paulis, phis)
      @test ok                                   # frame-dependent: k and natures (§4)
      worst = max(worst, max_deviation(ψ, dense_state(2, paulis, phis)))
      ψ, _ = run_truncating(T, 2, paulis, phis)
      worst = max(worst, max_deviation(ψ, dense_state(2, paulis, phis)))
    end
    @test worst <= EXACT_ATOL
  end

  @testset "angles outside (−π, π]" begin
    # Other multiples of π/4 and large angles: only the physical state is specified.
    worst = 0.0
    for ϕ in (5π/4, -5π/4, 7π/4, 2π, -2π, 9π/4, 3π, 3π + 0.3, -7.1, 12.0)
      paulis = [P"XY", P"ZX", P"-YI"]
      phis = [0.37, ϕ, 0.9]
      ψ, _, _, ok, _ = run_camps(T, 2, paulis, phis)
      # k bookkeeping (§4): the multiples of π/4 here are Clifford on QubitCAMPS (k unchanged);
      # on the old type `recognised_clifford` mirrors its |ϕ| ≤ π range.
      @test ok
      worst = max(worst, max_deviation(ψ, dense_state(2, paulis, phis)))
      ψ, _ = run_truncating(T, 2, paulis, phis)
      worst = max(worst, max_deviation(ψ, dense_state(2, paulis, phis)))
    end
    @test worst <= EXACT_ATOL
  end
end

@testset "apply!(ψ, P, ϕ, thl): the threshold truncates, $(stackname(T))" for T in STATE_TYPES
  # CONVENTIONS.md §6, truncated regime: the error must not grow as the threshold is
  # tightened and must vanish without truncation. No absolute bound is specified.
  bonddim = DisentangleCAMPS.bonddim
  N = 5
  paulis, phis = random_circuit(Xoshiro(1108), N, 40, [:general])
  ref = dense_state(N, paulis, phis)
  function run(thl)
    ψ = fresh_state(T, N)
    for (P, ϕ) in zip(paulis, phis)
      @test apply!(ψ, P, ϕ, thl) === nothing
    end
    return max_deviation(ψ, ref), 1 - norm(mps_part(ψ))^2, bonddim(ψ)
  end
  results = [run(thl) for thl in (1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 1e-4, 1e-15)]
  errs = getindex.(results, 1)
  lost = getindex.(results, 2)                   # weight discarded by the SVD cutoff
  @test errs[end] <= EXACT_ATOL
  @test errs[1] > 1e-4                           # the loosest threshold really truncates
  @test nonincreasing(errs)               # error does not grow as thl shrinks
  @test lost[1] > 1e-4 && abs(lost[end]) <= EXACT_ATOL
  @test nonincreasing(lost)               # … and neither does the discarded weight
  # The ladder is informative: several thresholds truncate, by different amounts.
  @test count(>(1e-6), errs) >= 3
  @test length(unique(round.(errs[errs .> 1e-6]; sigdigits = 3))) >= 3
  @test results[end][3] == 4                     # full bond dimension at N = 5 (frame-dependent)
end

"""
A Pauli sum Σ cᵢ Pᵢ (signed strings Pᵢ, real cᵢ) in the native sum type of T.
`QubitCAMPS`: a `QubitPauliSum` built with `push!`, which folds the sign of a string into its
coefficient (the contract `rotation_gate` relies on). Each key is a fresh copy, because a
`QubitPauliSum` keeps the key vector it is given.
"""
native_sum(::Type{CAMPS}, coeffs, ops) = PauliSum(coeffs, Stabilizer(ops))
function native_sum(::Type{mt.QubitCAMPS}, coeffs, ops)
  S = mt.QubitPauliSum(nqubits(ops[1]))
  for (c, P) in zip(coeffs, ops)
    push!(S, mt.QubitPauli(copy(P)), c)
  end
  return S
end
sum_expectation(ψ::CAMPS, S) = expectation(ψ, S)
sum_expectation(ψ::mt.QubitCAMPS, S) = mt.expectation(ψ, S)

@testset "expectation: signed strings and sums, $(stackname(T))" for T in STATE_TYPES
  # CONVENTIONS.md §5: O is a signed Pauli string or a Pauli sum.
  N = 3
  rng = Xoshiro(1109)
  paulis, phis = random_circuit(rng, N, 10, [:general, :general, π/4])
  ψ, _ = run_camps(T, N, paulis, phis)
  dense = dense_state(N, paulis, phis)
  ref(str) = DR.expect(dense, str)

  for str in ("XYZ", "ZIZ", "IYI", "YYX")
    @test real(state_expectation(ψ, pauliop(str))) ≈ ref(str) atol = EXACT_ATOL
    @test real(state_expectation(ψ, -pauliop(str))) ≈ -ref(str) atol = EXACT_ATOL
    @test abs(ref(str)) > 1e-3                   # not a vanishing value
  end

  # Sum with real coefficients, one of the strings carrying a minus sign.
  strs = ["XYZ", "ZIZ", "IYI", "III"]
  coeffs = [0.5, -1.25, 2.0, 0.75]
  ops = [pauliop(strs[1]), -pauliop(strs[2]), pauliop(strs[3]), pauliop(strs[4])]
  expected = 0.5 * ref("XYZ") + 1.25 * ref("ZIZ") + 2.0 * ref("IYI") + 0.75
  ev = sum_expectation(ψ, native_sum(T, coeffs, ops))
  @test real(ev) ≈ expected atol = EXACT_ATOL
  @test abs(imag(ev)) <= EXACT_ATOL
end

# == QubitCAMPS-only: Clifford-angle recognition, guards and pp wrappers (T9 gaps 2, 3, 7) ====
#
# CONVENTIONS.md §2: both `apply!` methods take `atol` and `rtol` (defaults 1e-10 and 0); every
# multiple of π/4 (any sign, any magnitude) within `isapprox(ϕ, mπ/4; atol, rtol)` is a
# Clifford angle, applied as a frame update with k unchanged. Near-Clifford angles are checked
# through k (and the MPS/frame split), not amplitudes, so the ~1e-9 floor of the MPO route
# plays no part. Most checks run on |0⟩ with P = X at N = 1, where X is `:disentanglable`
# (Q = X has an x bit beyond the k = 0 magic qubits): an angle NOT taken as Clifford gives
# k = 1 and an angle taken as Clifford gives k = 0.

import PauliPropagation as pp
using CampsPP.DisentangleCAMPS: evolve_deepcliffords, rotation_gate
import ITensorMPS

const QT = mt.QubitCAMPS

"k after one X rotation on a fresh one-qubit state (the four-argument `apply!`)."
kafter(ϕ; kw...) = apply!(fresh_state(QT, 1), 0, P"X", ϕ; kw...)

@testset "Clifford-angle tolerance defaults (QubitCAMPS)" begin
  # Default atol = 1e-10, rtol = 0.
  # δ = 1e-9 is 10x above the default atol: not Clifford (k = 1); δ = 1e-11 is 10x below: Clifford.
  @test kafter(π/4 + 1e-9) == 1
  @test kafter(π/4 + 1e-11) == 0
  @test kafter(π/4 - 1e-9) == 1
  @test kafter(π/4 - 1e-11) == 0
  # explicit atol: 1e-9 < 1e-8 is Clifford; 1e-9 > 1e-10 is not
  @test kafter(π/4 + 1e-9; atol = 1e-8) == 0
  @test kafter(π/4 + 1e-9; atol = 1e-10) == 1
  # relative: ϕ = π/4 (1 + 1e-9) is off by 7.85e-10 > atol = 1e-10 and rtol = 0, so k = 1;
  # with rtol = 1e-8 the tolerance is 1e-8 · 0.785 = 7.85e-9 ≥ 7.85e-10, so Clifford.
  @test kafter(π/4 * (1 + 1e-9)) == 1
  @test kafter(π/4 * (1 + 1e-9); rtol = 1e-8) == 0
  @test kafter(π/4 * (1 + 1e-9); rtol = 1e-12) == 1
  # near zero the relative tolerance is negligible (rtol · |ϕ| ~ 1e-17): only atol counts
  @test kafter(1e-9) == 1
  @test kafter(1e-11) == 0

  # Other Clifford angles, any sign and magnitude (m = round(ϕ/(π/4))), both sides of the cutoff.
  for m in (0, 1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, 8, -8, 9, -12)
    c = m * π / 4
    @test kafter(c) == 0
    @test kafter(c + 1e-11) == 0
    @test kafter(c - 1e-11) == 0
    @test kafter(c + 1e-9) == 1
    @test kafter(c - 1e-9) == 1
    @test kafter(c + 1e-9; atol = 1e-8) == 0
    @test kafter(c - 1e-9; atol = 1e-8) == 0
    # relative tolerance, meaningful away from 0 (|c|·1e-9 ≥ 7.85e-10 > atol for |c| ≥ π/4;
    # at |c| = 12·π/4 the relative offset is c·1e-9 ≈ 9e-9 < rtol·|c| = 9e-8)
    m == 0 && continue
    @test kafter(c * (1 + 1e-9)) == 1
    @test kafter(c * (1 + 1e-9); rtol = 1e-8) == 0
  end
end

@testset "Clifford-angle keyword forwarding (QubitCAMPS)" begin
  # Same near-Clifford one-gate circuit through every entry point: k = 1 by default,
  # k = 0 with the tolerance widened (atol = 1e-8, or rtol = 1e-8 for the relative case).
  cases = ((π/4 + 1e-9, (atol = 1e-8,)), (π/4 * (1 + 1e-9), (rtol = 1e-8,)))
  for (ϕ, kw) in cases
    @test evolve(fresh_state(QT, 1), 1, [P"X"], [ϕ])[2] == 1
    @test evolve(fresh_state(QT, 1), 1, [P"X"], [ϕ]; kw...)[2] == 0
    @test evolve_bonddim(fresh_state(QT, 1), 8, [P"X"], [ϕ])[2] == 1
    @test evolve_bonddim(fresh_state(QT, 1), 8, [P"X"], [ϕ]; kw...)[2] == 0
  end

  # evolve_deepcliffords: a uniformly random Clifford R precedes the gate. From the identity
  # frame, C ← R·C gives C = R, so Q = C†XX C = R†XX R is a random Pauli; with k = 0, k becomes
  # 1 exactly when Q has an x bit (§4, `spec_nature`), else stays 0. Widened tolerance: k = 0
  # always (Clifford gate). Default: k pinned per seed. R is regenerated by replaying the
  # seed: implementation-specific, per the docstring (layer s draws
  # `random_clifford(rng, N)` from the rng passed in), as in TestEntropies.jl. This catches a
  # k that is too high (e.g. every non-Clifford gate counted as a new magic qubit), which
  # leaves the state exact and so is invisible to the dense checks.
  N = 2
  kdefault, kexpected = Int[], Int[]
  for seed in 1:12, (ϕ, kw) in cases
    R = random_clifford(Xoshiro(seed), N)
    Q = (inv(R) * Stabilizer([P"XX"]))[1]       # R†·XX·R
    push!(kexpected, spec_nature(0, Q) == :disentanglable ? 1 : 0)
    push!(kdefault, evolve_deepcliffords(Xoshiro(seed), fresh_state(QT, N), 1, [P"XX"], [ϕ])[2])
    @test evolve_deepcliffords(Xoshiro(seed), fresh_state(QT, N), 1, [P"XX"], [ϕ]; kw...)[2] == 0
  end
  @test kdefault == kexpected
  # Not vacuous: both outcomes occur among the seeds.
  @test any(==(0), kexpected) && any(==(1), kexpected)

  # Truncating apply!(ψ, P, ϕ, thl; atol, rtol) returns nothing: observe it through the state.
  # Clifford path: the MPS is untouched (still |0⟩, overlap 1) and the frame changed (a
  # quarter turn about X conjugates Z → ∓Y, so the frame is not the identity).
  # Generic path: MPS becomes exp(iϕX)|0⟩, overlap with |0⟩ is cos ϕ ≈ 0.707, frame untouched.
  id1 = one(CliffordOperator, 1)
  for (ϕ, kw) in cases
    ψ = fresh_state(QT, 1)
    before = copy(mps_part(ψ))
    @test apply!(ψ, P"X", ϕ, 1e-15; kw...) === nothing
    @test abs(ITensorMPS.inner(before, mps_part(ψ))) ≈ 1 atol = 1e-12
    @test frame_dagger(ψ) != id1

    ψ = fresh_state(QT, 1)
    before = copy(mps_part(ψ))
    @test apply!(ψ, P"X", ϕ, 1e-15) === nothing                       # defaults: not Clifford
    @test abs(ITensorMPS.inner(before, mps_part(ψ))) ≈ cos(π/4) atol = 1e-6
    @test frame_dagger(ψ) == id1
  end
end

@testset "multiples of π/4 beyond ±π: k and frame (QubitCAMPS)" begin
  # exp(iϕX)|0⟩ for ϕ = mπ/4: k stays 0 (frame update) for every m, including |m| > 4.
  # Quarter turns n = m mod 4: n = 0 leaves the frame exactly the identity (a global phase
  # only); n ≠ 0 changes it; the MPS stays |0⟩ in all cases. State from the dense reference.
  id1 = one(CliffordOperator, 1)
  for m in (5, 6, 7, 8, 9, 12, 13, -5, -6, -7, -8, -12, 17, -17, 40, -40)
    ϕ = m * π / 4
    ψ = fresh_state(QT, 1)
    before = copy(mps_part(ψ))
    k = apply!(ψ, 0, P"X", ϕ)
    @test k == 0
    @test abs(ITensorMPS.inner(before, mps_part(ψ))) ≈ 1 atol = 1e-12
    @test (frame_dagger(ψ) == id1) == (mod(m, 4) == 0)
    @test max_deviation(ψ, dense_state(1, [P"X"], [ϕ])) <= EXACT_ATOL
  end
  # Same inside a two-qubit circuit (the k bookkeeping of `run_camps` is live there).
  for m in (5, 6, 7, 8, 9, 12, -5, -6, -7, -8, -12)
    paulis = [P"XY", P"ZX", P"-YI"]
    phis = [0.37, m * π / 4, 0.9]
    ψ, k, natures, ok, _ = run_camps(QT, 2, paulis, phis)
    @test ok
    @test natures[2] == :disentanglable   # the tested gate really is a Clifford-vs-magic decision
    @test max_deviation(ψ, dense_state(2, paulis, phis)) <= EXACT_ATOL
  end
end

@testset "guards: non-Hermitian strings and non-finite angles (QubitCAMPS)" begin
  # exp(iϕP) = cos ϕ + i sin ϕ P is unitary only for Hermitian P (§2): a phase ±i string
  # is rejected. The spec leaves the exception type open (spec_additions_needed).
  badP = (im * P"X", -im * P"XZ", im * P"ZZ")
  for P in badP, ϕ in (0.3, π/4)
    @test_throws ArgumentError apply!(fresh_state(QT, nqubits(P)), 0, P, ϕ)
    @test_throws ArgumentError apply!(fresh_state(QT, nqubits(P)), P, ϕ, 1e-15)
  end
  for P in badP
    @test_throws ArgumentError rotation_gate(P, 0.3)
  end
  @test rotation_gate(P"-XZ", 0.3) isa mt.QubitPauliSum     # ±1 signs are fine

  # Non-finite angles: rejected, with ArgumentError (not the InexactError of round(Int, NaN)).
  for ϕ in (NaN, Inf, -Inf)
    @test_throws ArgumentError apply!(fresh_state(QT, 1), 0, P"X", ϕ)
    @test_throws ArgumentError apply!(fresh_state(QT, 1), P"X", ϕ, 1e-15)
  end
end

@testset "evolve wrappers on pp gates vs dense (QubitCAMPS)" begin
  # `evolve` / `evolve_bonddim` taking vectors of pp.PauliString (signed) and pp.PauliRotation
  # must give the dense-reference state and the same k as the PauliOperator method.
  N = 3
  letters_of(P) = [x ? (z ? :Y : :X) : (z ? :Z : :I) for (x, z) in zip(xbit(P), zbit(P))]
  function ppstring(P)
    ls = letters_of(P)
    qs = [j for j in 1:N if ls[j] != :I]
    return pp.PauliString(N, ls[qs], qs, P.phase[] == 0x2 ? -1 : 1)
  end
  function pprot(P)
    ls = letters_of(P)
    qs = [j for j in 1:N if ls[j] != :I]
    return pp.PauliRotation(ls[qs], qs)
  end

  rng = Xoshiro(1121)
  paulis, phis = random_circuit(rng, N, 8, [:general, :general, π/4])
  # a single-qubit string on qubit N first, so a wrong qubit count in a wrapper is seen
  paulis[1] = P"IIX"
  strs = ppstring.(paulis)
  kref = run_camps(QT, N, paulis, phis)[2]
  dense = dense_state(N, paulis, phis)

  ψ, k = evolve(fresh_state(QT, N), length(paulis), strs, phis)
  @test k == kref
  @test max_deviation(ψ, dense) <= EXACT_ATOL
  ψ, k, s = evolve_bonddim(fresh_state(QT, N), 64, strs, phis)
  @test k == kref && s == length(paulis)
  @test max_deviation(ψ, dense) <= EXACT_ATOL

  # rotations carry no sign: use the unsigned circuit
  upaulis = [PauliOperator(0x0, xbit(P), zbit(P)) for P in paulis]
  rots = pprot.(upaulis)
  ukref = run_camps(QT, N, upaulis, phis)[2]
  udense = dense_state(N, upaulis, phis)
  ψ, k, s = evolve_bonddim(fresh_state(QT, N), 64, rots, phis)
  @test k == ukref && s == length(upaulis)
  @test max_deviation(ψ, udense) <= EXACT_ATOL
  ψ, k = evolve(fresh_state(QT, N), length(upaulis), rots, phis)
  @test k == ukref
  @test max_deviation(ψ, udense) <= EXACT_ATOL

  # prefix: evolve with t < length applies the first t gates only
  ψ, k = evolve(fresh_state(QT, N), 3, strs, phis)
  @test max_deviation(ψ, dense_state(N, paulis[1:3], phis[1:3])) <= EXACT_ATOL

  # keyword forwarding through every wrapper: near-Clifford single gate (k = 1 default, 0 widened)
  one_s = [pp.PauliString(1, [:X], [1], 1)]
  one_r = [pp.PauliRotation([:X], [1])]
  for (ϕ, kw) in ((π/4 + 1e-9, (atol = 1e-8,)), (π/4 * (1 + 1e-9), (rtol = 1e-8,)))
    for gates in (one_s, one_r)
      @test evolve(fresh_state(QT, 1), 1, gates, [ϕ])[2] == 1
      @test evolve(fresh_state(QT, 1), 1, gates, [ϕ]; kw...)[2] == 0
      @test evolve_bonddim(fresh_state(QT, 1), 8, gates, [ϕ])[2] == 1
      @test evolve_bonddim(fresh_state(QT, 1), 8, gates, [ϕ]; kw...)[2] == 0
    end
  end
end
