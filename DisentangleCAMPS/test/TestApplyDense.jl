# Randomized-exact tests of `apply!` and `evolve*` against the dense reference (T4/T6).
# Seeded random circuits; each instance is compared exactly, through all 4^N Pauli
# expectation values, so the comparison is frame-invariant (CONVENTIONS.md §4).
#
# Expected values come from DenseReference.jl only. The code's own `paulinature` is used
# for one thing: counting which branches the random circuits reached.

using DisentangleCAMPS
using DisentangleCAMPS: paulinature, evolve_bonddim
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
function max_deviation(ψ::CAMPS, dense)
  N = length(ψ)
  dev = 0.0
  for (str, ref) in zip(DR.all_paulis(N), DR.all_expectations(dense))
    ev = expectation(ψ, pauliop(str))
    dev = max(dev, abs(real(ev) - ref), abs(imag(ev)))
  end
  return dev
end

dense_state(N, paulis, phis) = DR.run_circuit(DR.zero_state(N), paulistr.(paulis), phis)

# Angles the old code recognises as Clifford (CONVENTIONS.md §2); such gates leave k unchanged.
const OLD_CLIFFORD_ANGLES = (0.0, π/4, π/2, 3π/4, π)
recognised_clifford(ϕ) = any(isapprox(abs(ϕ), c; atol = 1e-9) for c in OLD_CLIFFORD_ANGLES)

"""
Run the 4-argument `apply!` over the circuit from |0…0⟩.
Returns (ψ, k, natures, bookkeeping_ok, n_trivial_with_k>0).

Frame-dependent bookkeeping, checked gate by gate against CONVENTIONS.md §4 with the
independent `spec_nature`: Q = C†PC = Cdag·P·Cdag⁻¹ (old frame), `paulinature` must agree
with the spec rule, and k must grow by one exactly for a `:disentanglable` gate at an
angle not recognised as Clifford.
"""
function run_camps(N, paulis, phis)
  ψ = CAMPS(N)
  k = 0
  kspec = 0
  natures = Symbol[]
  ok = true
  ntriv = 0
  for (P, ϕ) in zip(paulis, phis)
    nature = spec_nature(k, ψ.Cdag * P)
    push!(natures, nature)
    ntriv += nature == :trivial && k > 0
    ok &= paulinature(k, inv(ψ.Cdag), P) == nature
    kspec += nature == :disentanglable && !recognised_clifford(ϕ)
    k = apply!(ψ, k, P, ϕ)
    ok &= k == kspec
  end
  return ψ, k, natures, ok, ntriv
end

"Seed for a named test set: stable across Julia versions, different for every name."
setseed(base, name) = Xoshiro(base + sum(Int, codeunits(name)))

"Worst deviation over `ncirc` random circuits, plus the fraction of circuits within tolerance."
function sweep(rng, N, depth, ncirc, angle_sources; runner = run_camps)
  worst = 0.0
  nok = 0
  natures = Symbol[]
  bookkeeping = true
  ntriv = 0
  for _ in 1:ncirc
    paulis, phis = random_circuit(rng, N, depth, angle_sources)
    ψ, k, nat, ok, nt = runner(N, paulis, phis)
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

@testset "apply!(ψ, k, P, ϕ) vs dense" begin
  @testset "general angles" begin
    rng = Xoshiro(1101)
    natures = Symbol[]
    ntrivial_kpos = 0
    for (N, depth, ncirc) in [(2, 8, 40), (3, 10, 60), (4, 12, 10)]
      r = sweep(rng, N, depth, ncirc, [:general])
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
    ψ, k, _, ok, _ = run_camps(3, paulis, phis)
    @test k == 3
    @test ok
    @test max_deviation(ψ, dense_state(3, paulis, phis)) <= EXACT_ATOL
  end

  # Clifford branches that agree with exp(iϕP), alone and mixed with general angles.
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
      r = sweep(rng, N, depth, ncirc, sources)
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
      @test r.bookkeeping                        # frame-dependent: natures and k
    end
  end

  # Known bugs in the old `apply!` (plan, hard spot 6). Each line must turn green after
  # Phase 2: replace `@test_broken` by `@test` then.
  @testset "known special-angle bugs: $name" for (name, sources) in [
    ("−π/4 only", [-π/4]),
    ("+3π/4 only", [3π/4]),
    ("±π/4 mixed", [π/4, -π/4]),
    ("general with −π/4", [:general, -π/4]),
    ("general with +3π/4", [:general, 3π/4]),
    ("general with ±π/2", [:general, π/2, -π/2]),
    ("+π/4 with ±π/2", [π/4, π/2, -π/2]),
  ]
    rng = setseed(1104, name)
    r = sweep(rng, 3, 10, 20, sources)
    @test_broken r.nok == r.n
  end
end

@testset "apply!(ψ, P, ϕ, thl) vs dense (no truncation)" begin
  # The truncating method with a negligible cutoff must be exact too.
  function run_truncating(N, paulis, phis)
    ψ = CAMPS(N)
    for (P, ϕ) in zip(paulis, phis)
      apply!(ψ, P, ϕ, 1e-15)
    end
    return ψ, 0, Symbol[], true, 0
  end

  @testset "$name" for (name, sources) in [
    ("general", [:general]),
    ("general with +π/4, −3π/4, 0, π", [:general, :general, π/4, -3π/4, 0.0, π]),
    # Without π/4 gates this method never changes the frame except by Paulis, and the
    # ±π/2 branch is then correct (Paulis commute as Clifford operators).
    ("general with ±π/2", [:general, π/2, -π/2]),
  ]
    rng = setseed(1105, name)
    for (N, depth, ncirc) in [(2, 8, 20), (3, 10, 20), (4, 12, 5)]
      r = sweep(rng, N, depth, ncirc, sources; runner = run_truncating)
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
    end
  end

  @testset "known special-angle bugs: $name" for (name, sources) in [
    ("general with −π/4", [:general, -π/4]),
    ("general with +3π/4", [:general, 3π/4]),
    ("+π/4 with ±π/2", [π/4, π/2, -π/2]),
    ("general with +π/4 and ±π/2", [:general, π/4, π/2, -π/2]),
  ]
    rng = setseed(1106, name)
    r = sweep(rng, 3, 10, 20, sources; runner = run_truncating)
    @test_broken r.nok == r.n
  end
end

@testset "evolve and evolve_bonddim vs dense" begin
  bonddim = DisentangleCAMPS.bonddim
  rng = Xoshiro(1107)
  N = 3
  paulis, phis = random_circuit(rng, N, 12, [:general, :general, π/4])

  # evolve(ψ, t, …) applies the first t gates and returns the same k as gate-by-gate `apply!`.
  for t in (0, 5, 12)
    ψ, k = evolve(CAMPS(N), t, paulis, phis)
    @test max_deviation(ψ, dense_state(N, paulis[1:t], phis[1:t])) <= EXACT_ATOL
    _, kref, _, ok, _ = run_camps(N, paulis[1:t], phis[1:t])
    @test ok && k == kref
  end

  # evolve_bonddim with χ above the maximum possible bond dimension (2 for N = 3) runs
  # the whole circuit and reports s = number of gates applied.
  ψ, k, s = evolve_bonddim(CAMPS(N), 64, paulis, phis)
  @test s == length(paulis)
  @test k == run_camps(N, paulis, phis)[2]
  @test max_deviation(ψ, dense_state(N, paulis, phis)) <= EXACT_ATOL

  # With χ = 2 it stops at the first gate after which the bond dimension is ≥ χ
  # (frame-dependent): reached after s gates, not yet after s − 1.
  ψ, k, s = evolve_bonddim(CAMPS(N), 2, paulis, phis)
  @test 1 < s < length(paulis)
  @test bonddim(ψ) >= 2
  @test bonddim(evolve(CAMPS(N), s - 1, paulis, phis)[1]) < 2
  @test k == run_camps(N, paulis[1:s], phis[1:s])[2]
  @test max_deviation(ψ, dense_state(N, paulis[1:s], phis[1:s])) <= EXACT_ATOL
end

@testset "evolve_bonddim: hand-built stopping point" begin
  # N = 2 from |00⟩. Gates 1 and 2 are single-qubit X rotations: each is `:disentanglable`
  # with trivial disentangler (no other qubit involved), so the frame stays the identity,
  # the MPS is the product exp(i·0.4·X)|0⟩ ⊗ exp(i·0.7·X)|0⟩ (bond dimension 1) and k = 2.
  # Gate 3, exp(i·0.5·ZZ), is `:logical` and entangles the two qubits: bond dimension 2.
  # (ZZ is not a product operator on that state: ⟨ZZ⟩ = cos 0.8 · cos 1.4 ≠ ±1.)
  bonddim = DisentangleCAMPS.bonddim
  paulis = [P"XI", P"IX", P"ZZ", P"XY"]
  phis = [0.4, 0.7, 0.5, 0.3]

  ψ, k, s = evolve_bonddim(CAMPS(2), 2, paulis, phis)
  @test s == 3                                   # stops right after the entangling gate
  @test k == 2
  @test bonddim(ψ) == 2
  @test max_deviation(ψ, dense_state(2, paulis[1:3], phis[1:3])) <= EXACT_ATOL

  ψ, k, s = evolve_bonddim(CAMPS(2), 3, paulis, phis)   # χ = 3 cannot be reached at N = 2
  @test s == 4
  @test k == 2
  @test max_deviation(ψ, dense_state(2, paulis, phis)) <= EXACT_ATOL

  ψ, k, s = evolve_bonddim(CAMPS(2), 1, paulis, phis)   # already at χ: nothing is applied
  @test s == 0 && k == 0
end

@testset "angles next to Clifford angles" begin
  # CONVENTIONS.md §2: the state must be exp(iϕP)|ψ⟩ for every real ϕ. An implementation
  # that rounds ϕ = c ± 1e-3 to the Clifford angle c is wrong by about 2e-3.
  multiples = [m * π / 4 for m in -4:4]

  @testset "N = 1: state and k" begin
    for c in multiples, δ in (-1e-3, 1e-3)
      ϕ = c + δ
      ψ = CAMPS(1)
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
      ψ, _ = run_camps(2, paulis, phis)
      worst = max(worst, max_deviation(ψ, dense_state(2, paulis, phis)))
      ψ = CAMPS(2)
      for (P, ϕ) in zip(paulis, phis)
        apply!(ψ, P, ϕ, 1e-15)                   # truncating method
      end
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
      ψ, _ = run_camps(2, paulis, phis)
      worst = max(worst, max_deviation(ψ, dense_state(2, paulis, phis)))
    end
    @test worst <= EXACT_ATOL
  end
end

@testset "apply!(ψ, P, ϕ, thl): the threshold truncates" begin
  # CONVENTIONS.md §6, truncated regime: the error must not grow as the threshold is
  # tightened and must vanish without truncation. No absolute bound is specified.
  bonddim = DisentangleCAMPS.bonddim
  N = 5
  paulis, phis = random_circuit(Xoshiro(1108), N, 40, [:general])
  ref = dense_state(N, paulis, phis)
  function run(thl)
    ψ = CAMPS(N)
    for (P, ϕ) in zip(paulis, phis)
      @test apply!(ψ, P, ϕ, thl) === nothing
    end
    return max_deviation(ψ, ref), 1 - norm(ψ)^2, bonddim(ψ)
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

@testset "expectation: signed strings and sums" begin
  # CONVENTIONS.md §5: O is a signed Pauli string or a Pauli sum.
  N = 3
  rng = Xoshiro(1109)
  paulis, phis = random_circuit(rng, N, 10, [:general, :general, π/4])
  ψ, _ = run_camps(N, paulis, phis)
  dense = dense_state(N, paulis, phis)
  ref(str) = DR.expect(dense, str)

  for str in ("XYZ", "ZIZ", "IYI", "YYX")
    @test real(expectation(ψ, pauliop(str))) ≈ ref(str) atol = EXACT_ATOL
    @test real(expectation(ψ, -pauliop(str))) ≈ -ref(str) atol = EXACT_ATOL
    @test abs(ref(str)) > 1e-3                   # not a vanishing value
  end

  # Sum with real coefficients, one of the strings carrying a minus sign.
  strs = ["XYZ", "ZIZ", "IYI", "III"]
  coeffs = [0.5, -1.25, 2.0, 0.75]
  ops = Stabilizer([pauliop(strs[1]), -pauliop(strs[2]), pauliop(strs[3]), pauliop(strs[4])])
  expected = 0.5 * ref("XYZ") + 1.25 * ref("ZIZ") + 2.0 * ref("IYI") + 0.75
  ev = expectation(ψ, PauliSum(coeffs, ops))
  @test real(ev) ≈ expected atol = EXACT_ATOL
  @test abs(imag(ev)) <= EXACT_ATOL
end
