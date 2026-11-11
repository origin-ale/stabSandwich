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

"Run the 4-argument `apply!` over the circuit from |0…0⟩. Returns (ψ, k, natures seen)."
function run_camps(N, paulis, phis)
  ψ = CAMPS(N)
  k = 0
  natures = Symbol[]
  for (P, ϕ) in zip(paulis, phis)
    push!(natures, paulinature(k, inv(ψ.Cdag), P))
    k = apply!(ψ, k, P, ϕ)
  end
  return ψ, k, natures
end

"Worst deviation over `ncirc` random circuits, plus the fraction of circuits within tolerance."
function sweep(rng, N, depth, ncirc, angle_sources; runner = run_camps)
  worst = 0.0
  nok = 0
  natures = Symbol[]
  for _ in 1:ncirc
    paulis, phis = random_circuit(rng, N, depth, angle_sources)
    ψ, k, nat = runner(N, paulis, phis)
    append!(natures, nat)
    dev = max_deviation(ψ, dense_state(N, paulis, phis))
    worst = max(worst, dev)
    nok += dev <= EXACT_ATOL
    @assert 0 <= k <= N
  end
  return (worst = worst, nok = nok, n = ncirc, natures = natures)
end

# -- Tests ------------------------------------------------------------------------------

@testset "apply!(ψ, k, P, ϕ) vs dense" begin
  @testset "general angles" begin
    rng = Xoshiro(1101)
    natures = Symbol[]
    for (N, depth, ncirc) in [(2, 8, 30), (3, 10, 30), (4, 12, 8)]
      r = sweep(rng, N, depth, ncirc, [:general])
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
      append!(natures, r.natures)
    end
    # Every `paulinature` branch was reached (no vacuous coverage).
    @test count(==(:disentanglable), natures) >= 50
    @test count(==(:logical), natures) >= 50
    @test count(==(:trivial), natures) >= 5
  end

  @testset "k reaches N" begin
    # After enough general rotations every qubit is a magic qubit.
    rng = Xoshiro(1102)
    paulis, phis = random_circuit(rng, 3, 12, [:general])
    ψ, k, _ = run_camps(3, paulis, phis)
    @test k == 3
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
    rng = Xoshiro(1103)
    for (N, depth, ncirc) in [(2, 8, 20), (3, 10, 20)]
      r = sweep(rng, N, depth, ncirc, sources)
      @test r.worst <= EXACT_ATOL
      @test r.nok == r.n
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
    rng = Xoshiro(1104)
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
    return ψ, 0, Symbol[]
  end

  @testset "$name" for (name, sources) in [
    ("general", [:general]),
    ("general with +π/4, −3π/4, 0, π", [:general, :general, π/4, -3π/4, 0.0, π]),
    # Without π/4 gates this method never changes the frame except by Paulis, and the
    # ±π/2 branch is then correct (Paulis commute as Clifford operators).
    ("general with ±π/2", [:general, π/2, -π/2]),
  ]
    rng = Xoshiro(1105)
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
    rng = Xoshiro(1106)
    r = sweep(rng, 3, 10, 20, sources; runner = run_truncating)
    @test_broken r.nok == r.n
  end
end

@testset "evolve and evolve_bonddim vs dense" begin
  rng = Xoshiro(1107)
  N = 3
  paulis, phis = random_circuit(rng, N, 12, [:general, :general, π/4])

  # evolve(ψ, t, …) applies the first t gates.
  for t in (0, 5, 12)
    ψ, k = evolve(CAMPS(N), t, paulis, phis)
    @test max_deviation(ψ, dense_state(N, paulis[1:t], phis[1:t])) <= EXACT_ATOL
  end

  # evolve_bonddim with χ above the maximum possible bond dimension (2 for N = 3) runs
  # the whole circuit and reports s = number of gates applied.
  ψ, k, s = evolve_bonddim(CAMPS(N), 64, paulis, phis)
  @test s == length(paulis)
  @test max_deviation(ψ, dense_state(N, paulis, phis)) <= EXACT_ATOL

  # With χ = 2 it stops as soon as the bond dimension reaches 2; the state must equal
  # the dense state after exactly s gates.
  ψ, k, s = evolve_bonddim(CAMPS(N), 2, paulis, phis)
  @test 0 < s <= length(paulis)
  @test max_deviation(ψ, dense_state(N, paulis[1:s], phis[1:s])) <= EXACT_ATOL
end
