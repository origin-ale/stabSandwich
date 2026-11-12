# Clifford-angle circuits at N = 46: `apply!` against the stabilizer reference (T6, large N).
# Random signed Pauli strings (about 3N/4 non-identity letters) at angles that are multiples
# of π/4. Expected values come from StabilizerReference.jl (CONVENTIONS.md §2), which is
# itself checked against the dense reference below.
#
# Frame-invariant comparison: ⟨g⟩ = 1 for all N generators g of the reference state, and
# ‖ψ‖² = 1. Together these fix the state exactly. Two guards on `expectation` itself:
# random-signed products of generators must give exactly their sign, ±1 (an `expectation`
# that drops signs or takes abs would otherwise pass, and with it the ±π/4 sign bugs, which
# only turn g into −g); random strings outside the stabilizer group must give 0.
#
# Reuses `pauliop`, `paulistr`, `random_gate_pauli`, `random_circuit`, `dense_state` and
# `setseed` from TestApplyDense.jl, which must be included before this file.

isdefined(@__MODULE__, :random_circuit) ||
  error("TestCliffordLargeN.jl needs the helpers of TestApplyDense.jl: include that file first")

isdefined(@__MODULE__, :StabilizerReference) || include("StabilizerReference.jl")
const SR = StabilizerReference

const N_LARGE = 46
const PI4_MULTIPLES = [k * π / 4 for k in -4:4]

@testset "stabilizer reference vs dense (self-test)" begin
  rng = Xoshiro(1201)
  for (N, ones) in [(2, Int[]), (3, Int[]), (3, [1, 3])]
    for _ in 1:15
      paulis, phis = random_circuit(rng, N, 10, PI4_MULTIPLES)
      gens = SR.run_generators(SR.basis_generators(N, ones), paulis, phis)
      dense = DR.run_circuit(DR.basis_state(N, ones), paulistr.(paulis), phis)
      # Every generator has ⟨g⟩ = 1 on the dense state …
      @test all(abs(DR.expect(dense, paulistr(g)) - 1) <= EXACT_ATOL for g in gens)
      # … and the reference's `expect` agrees with dense on all 4^N strings.
      S = Stabilizer(gens)
      @test maximum(abs(SR.stab_expect(S, pauliop(str)) - ref)
                    for (str, ref) in zip(DR.all_paulis(N), DR.all_expectations(dense))) <= EXACT_ATOL
    end
  end
  # Sign-sensitive single cases (CONVENTIONS.md §2): exp(iπ/4·X)|0⟩ = (|0⟩ + i|1⟩)/√2 has
  # ⟨Y⟩ = +1, since Z ↦ iXZ = i(−iY) = Y …
  @test SR.run_generators(SR.basis_generators(1, Int[]), [P"X"], [π/4]) == [P"Y"]
  # … the opposite angle gives −Y, ±π/2 flips the sign of Z, and 0, π change nothing.
  @test SR.run_generators(SR.basis_generators(1, Int[]), [P"X"], [-π/4]) == [P"-Y"]
  @test SR.run_generators(SR.basis_generators(1, Int[]), [P"X"], [π/2]) == [P"-Z"]
  @test SR.run_generators(SR.basis_generators(1, Int[]), [P"X"], [π]) == [P"Z"]
  @test_throws ArgumentError SR.rotate_generators!(SR.basis_generators(1, Int[]), P"X", 0.3)

  # At N = 46 (QuantumClifford's multi-word bit packing): the inverse circuit returns the
  # generators to the start, and they stay pairwise commuting and independent.
  ok = true
  for _ in 1:5
    paulis, phis = random_circuit(rng, N_LARGE, 150, PI4_MULTIPLES)
    g0 = SR.basis_generators(N_LARGE, findall(rand(rng, Bool, N_LARGE)))
    g = SR.run_generators(g0, paulis, phis)
    ok &= SR.run_generators(g, reverse(paulis), -reverse(phis)) == g0
    ok &= g != g0
    ok &= !any(SR.anticommutes(a, b) for a in g, b in g)
    ok &= QuantumClifford.rank(MixedDestabilizer(Stabilizer(g))) == N_LARGE
  end
  @test ok
end

"""
Largest deviation of the CAMPS ψ from the stabilizer state with generators `gens`:
|⟨g⟩ − 1| (and |Im⟨g⟩|) over all generators, |‖ψ‖² − 1|, |⟨σ·G⟩ − σ| over `nprobe`
products G of random subsets of generators with a random sign σ (⟨G⟩ = 1 because the
generators commute and each is +1), and |⟨A⟩| over `nprobe` random strings A that
anticommute with some generator (so ⟨A⟩ = 0 on the reference state).
"""
function stab_deviation(ψ::CAMPS, gens, rng; nprobe = 10)
  N = length(gens)
  dev = abs(norm(ψ.mps)^2 - 1)
  for g in gens
    ev = expectation(ψ, g)
    dev = max(dev, abs(real(ev) - 1), abs(imag(ev)))
  end
  for _ in 1:nprobe
    G = prod(gens[randsubseq(rng, 1:N, 0.5)]; init = PauliOperator(0x0, falses(N), falses(N)))
    G.phase[] in (0x0, 0x2) || error("product of commuting generators is not Hermitian")
    σ = rand(rng, (-1, 1))
    ev = expectation(ψ, σ == 1 ? G : -G)
    dev = max(dev, abs(real(ev) - σ), abs(imag(ev)))
  end
  nprobed = 0
  while nprobed < nprobe
    A = random_gate_pauli(rng, N)
    any(SR.anticommutes(A, g) for g in gens) || continue
    dev = max(dev, abs(expectation(ψ, A)))
    nprobed += 1
  end
  return dev
end

meanweight(gens) = sum(count(xbit(g) .| zbit(g)) for g in gens) / length(gens)

"""
`ncirc` random circuits of `depth` gates at N = 46 from |0…0⟩, run with the 4-argument
`apply!` (or the truncating one). Returns the worst deviation, the number of circuits
within tolerance, whether k and the bond dimension stayed trivial (frame-dependent: every
gate is at a Clifford angle), and the smallest mean generator weight reached.
"""
function large_sweep(rng, depth, ncirc, sources; truncating = false)
  worst, nok, trivial, minweight = 0.0, 0, true, Inf
  for _ in 1:ncirc
    paulis, phis = random_circuit(rng, N_LARGE, depth, sources)
    ψ = CAMPS(N_LARGE)
    k = 0
    for (P, ϕ) in zip(paulis, phis)
      truncating ? apply!(ψ, P, ϕ, 1e-15) : (k = apply!(ψ, k, P, ϕ))
    end
    trivial &= k == 0 && DisentangleCAMPS.bonddim(ψ) == 1
    gens = SR.run_generators(SR.basis_generators(N_LARGE, Int[]), paulis, phis)
    minweight = min(minweight, meanweight(gens))
    dev = stab_deviation(ψ, gens, rng)
    worst = max(worst, dev)
    nok += dev <= EXACT_ATOL
  end
  return (; worst, nok, n = ncirc, trivial, minweight)
end

@testset "apply! at N = 46, Clifford angles, vs stabilizer reference" begin
  @testset "$(truncating ? "apply!(ψ, P, ϕ, thl)" : "apply!(ψ, k, P, ϕ)"): $name" for
      truncating in (false, true),
      (name, sources) in [
        ("+π/4, −3π/4, 0, ±π", [π/4, π/4, -3π/4, 0.0, π, -π]),
        ("±π/2, 0, π", [π/2, -π/2, 0.0, π]),   # Pauli frame only: the ±π/2 branch is correct
      ]
    rng = setseed(1202 + truncating, name)
    r = large_sweep(rng, 120, 3, sources; truncating)
    @test r.worst <= EXACT_ATOL
    @test r.nok == r.n
    @test r.trivial                                # frame-dependent: k = 0, bond dimension 1
    # Not vacuous: π/4 gates spread the generators to high weight; π/2 gates only flip signs.
    occursin("π/4", name) ? (@test r.minweight > N_LARGE / 3) : (@test r.minweight == 1)
  end

  # Known bugs in the old `apply!` (CONVENTIONS.md §2, plan hard spot 6), at N = 46. They
  # must turn green after Phase 2: replace `@test_broken` by `@test` then.
  @testset "known special-angle bugs ($(truncating ? "truncating" : "4-argument")): $name" for
      truncating in (false, true),
      (name, sources) in [
        ("−π/4 only", [-π/4]),
        ("+3π/4 only", [3π/4]),
        ("+π/4 with ±π/2", [π/4, π/2, -π/2]),
      ]
    rng = setseed(1204 + truncating, name)
    r = large_sweep(rng, 60, 2, sources; truncating)
    @test_broken r.nok == r.n
  end
end
