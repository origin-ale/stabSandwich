# Entropy functions and `evolve_deepcliffords` (T6). Expected values: CONVENTIONS.md §5,
# computed from dense-reference states by the formulas written here, or by hand.

using DisentangleCAMPS
using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using LinearAlgebra: svdvals
using Random
using Test

isdefined(@__MODULE__, :DenseReference) || include("DenseReference.jl")

# -- Reference formulas (CONVENTIONS.md §5) ----------------------------------------------

"Dense state → MPS. The dense index has qubit 1 as its least significant bit, i.e. the first array index."
dense_to_mps(ψ, N) = MPS(reshape(ψ, ntuple(_ -> 2, N)), siteinds("Qubit", N); cutoff = 1e-15)

"von Neumann entropy, base 2, of the bipartition (1…b | b+1…N) of a normalised dense state."
function dense_bond_entropy(ψ, N, b)
  p = svdvals(reshape(ψ, 2^b, 2^(N - b))) .^ 2
  return -sum(x -> x > 0 ? x * log2(x) : 0.0, p)
end

"Stabilizer Rényi entropy M_α of a normalised dense state; α = 1 and the linear entropy (α = −1 in the code)."
function dense_sre(ψ, N, α)
  Ξ = DenseReference.all_expectations(ψ) .^ 2 ./ 2^N
  α == 1 && return -sum(x -> x > 0 ? x * log2(x) : 0.0, Ξ) - N
  α == -1 && return 1 - 2^N * sum(Ξ .^ 2)
  return log2(sum(Ξ .^ α)) / (1 - α) - N
end

h2(p) = -p * log2(p) - (1 - p) * log2(1 - p)     # binary entropy, base 2

const ENT_ATOL = 1e-10

@testset "reference formulas on hand cases" begin
  DRm = DenseReference
  # Single-qubit magic state exp(iπ/8·X)|0⟩: Bloch (0, 1/√2, 1/√2), Ξ = (1, 0, ½, ½)/2,
  # ΣΞ² = 3/8, so M₂ = −log₂(3/8) − 1 = log₂(4/3).
  m = DRm.rotate(DRm.zero_state(1), "X", π/8)
  @test dense_sre(m, 1, 2) ≈ log2(4/3) atol = 1e-12
  @test dense_sre(m, 1, -1) ≈ 1 - 2 * 3/8 atol = 1e-12          # linear entropy = 1/4
  # Bell state (|00⟩ + |11⟩)/√2: a stabilizer state with one ebit.
  bell = ComplexF64[1, 0, 0, 1] ./ sqrt(2)
  @test dense_sre(bell, 2, 2) ≈ 0 atol = 1e-12
  @test dense_sre(bell, 2, 1) ≈ 0 atol = 1e-12
  @test dense_bond_entropy(bell, 2, 1) ≈ 1 atol = 1e-12
end

@testset "sEntropy and eEntropys! on hand cases" begin
  DRm = DenseReference
  # |00⟩: both vanish.
  ψ = dense_to_mps(DRm.zero_state(2), 2)
  @test sEntropy(ψ; α = 2.0) ≈ 0 atol = ENT_ATOL
  @test eEntropys!(copy(ψ)) ≈ [0.0] atol = ENT_ATOL

  # Bell state: M_α = 0 for every α, entanglement entropy 1 (log base 2).
  bell = dense_to_mps(ComplexF64[1, 0, 0, 1] ./ sqrt(2), 2)
  @test sEntropy(bell; α = 2.0) ≈ 0 atol = ENT_ATOL
  @test sEntropy(bell; α = 1.0) ≈ 0 atol = ENT_ATOL
  @test sEntropy(bell; α = -1.0) ≈ 0 atol = ENT_ATOL
  @test eEntropys!(copy(bell)) ≈ [1.0] atol = ENT_ATOL

  # Magic qubit next to |0⟩: additive, so M₂ = log₂(4/3) + 0, not divided by N; no entanglement.
  magic = dense_to_mps(DRm.rotate(DRm.zero_state(2), "XI", π/8), 2)
  @test sEntropy(magic; α = 2.0) ≈ log2(4/3) atol = ENT_ATOL
  @test sEntropy(magic; α = -1.0) ≈ 0.25 atol = ENT_ATOL
  @test eEntropys!(copy(magic)) ≈ [0.0] atol = ENT_ATOL

  # cos θ|00⟩ + i sin θ|11⟩: Schmidt probabilities (cos²θ, sin²θ).
  θ = 0.3
  ent = dense_to_mps(DRm.rotate(DRm.zero_state(2), "XX", θ), 2)
  @test eEntropys!(copy(ent)) ≈ [h2(cos(θ)^2)] atol = ENT_ATOL

  # GHZ on 3 qubits: a stabilizer state, one ebit across each of the N − 1 bonds.
  ghz = zeros(ComplexF64, 8); ghz[1] = ghz[8] = 1 / sqrt(2)
  ghz_mps = dense_to_mps(ghz, 3)
  @test sEntropy(ghz_mps; α = 2.0) ≈ 0 atol = ENT_ATOL
  @test eEntropys!(copy(ghz_mps)) ≈ [1.0, 1.0] atol = ENT_ATOL
end

@testset "sEntropy and eEntropys! vs dense, random states" begin
  DRm = DenseReference
  rng = Xoshiro(1201)
  letters = ['I', 'X', 'Y', 'Z']
  for N in (2, 3, 4), _ in 1:4
    strs = [String(rand(rng, letters, N)) for _ in 1:3N]
    phis = [(2rand(rng) - 1) * π for _ in strs]
    dense = DRm.run_circuit(DRm.zero_state(N), strs, phis)
    ψ = dense_to_mps(dense, N)

    for α in (2.0, 1.0, -1.0, 3.0)
      @test sEntropy(ψ; α) ≈ dense_sre(dense, N, α) atol = 1e-9
    end
    @test sEntropy(ψ, [2.0, 1.0]) ≈ [dense_sre(dense, N, 2), dense_sre(dense, N, 1)] atol = 1e-9

    ee = eEntropys!(copy(ψ))
    @test length(ee) == N - 1                    # element b: bond between sites b and b+1
    @test ee ≈ [dense_bond_entropy(dense, N, b) for b in 1:N-1] atol = 1e-9
  end
end

@testset "sampled sEntropy agrees with the exact value" begin
  # The sampler draws from the global RNG, so it is seeded globally here (Julia 1.13).
  # Statistical tolerance: fixed seed, 4000 samples, tolerance 0.05 on values of order 0.5.
  DRm = DenseReference
  N = 3
  dense = DRm.run_circuit(DRm.zero_state(N), ["XYI", "IZX", "YIY", "ZXZ"], [0.4, 0.9, -0.6, 1.2])
  ψ = dense_to_mps(dense, N)
  exact = dense_sre(dense, N, 2)
  @test exact > 0.2                              # a genuinely magic state
  Random.seed!(1202)
  @test sEntropy(ψ, 4000; α = 2.0) ≈ exact atol = 0.05
end

@testset "stabilizer Rényi entropy is frame-invariant on a CAMPS" begin
  # CONVENTIONS.md §5: M_α of the MPS part equals M_α of the physical state.
  DRm = DenseReference
  N = 3
  rng = Xoshiro(1203)
  paulis, phis = random_circuit(rng, N, 10, [:general, :general, π/4])
  ψ, k, _ = run_camps(N, paulis, phis)
  @test ψ.Cdag != one(CliffordOperator, N)       # the frame is not trivial
  dense = dense_state(N, paulis, phis)
  @test sEntropy(ψ; α = 2.0) ≈ dense_sre(dense, N, 2) atol = 1e-9
  @test sEntropy(ψ.mps; α = 2.0) ≈ dense_sre(dense, N, 2) atol = 1e-9
end

# Implementation-specific: this test regenerates the random Cliffords by replaying the
# global RNG through the old generator, `random_clifford_circuit(N, 2N^2)`. A correct port
# that draws its Cliffords differently fails here by construction; it should take an `rng`
# (or the Cliffords) as an argument, and this test should then pass them in.
@testset "evolve_deepcliffords vs dense" begin
  # Layer s applies a random Clifford C_s (depth 2N², from the *global* RNG) and then the
  # rotation (P_s, ϕ_s). The Cliffords are not returned, so they are regenerated here by
  # replaying the global RNG with the same seed and the same generator call.
  #
  # Expected values. With G_s = C_s ⋯ C_1,
  #   |ψ⟩ = R_t C_t ⋯ R_1 C_1 |0⟩ = G_t · Π_s exp(iϕ_s · G_s† P_s G_s) |0⟩ =: G_t |φ⟩,
  # so ⟨ψ|P|ψ⟩ = ⟨φ| G_t† P G_t |φ⟩. |φ⟩ comes from the dense reference, and the
  # conjugations from QuantumClifford (`inv(G) * P` = G†PG, pinned in TestGoldenSet.jl).
  N, t = 3, 5
  rng = Xoshiro(1204)
  paulis, phis = random_circuit(rng, N, t, [:general])

  seed = 4321
  Random.seed!(seed)
  cliffords = [random_clifford_circuit(N, 2N^2) for _ in 1:t]
  Random.seed!(seed)
  ψ, k = evolve_deepcliffords(CAMPS(N), t, paulis, phis)

  G = one(CliffordOperator, N)
  rotated = PauliOperator[]
  for s in 1:t
    G = cliffords[s] * G                         # C_s acts after the earlier ones
    push!(rotated, inv(G) * paulis[s])
  end
  φ = DenseReference.run_circuit(DenseReference.zero_state(N), paulistr.(rotated), phis)

  worst = 0.0
  for str in DenseReference.all_paulis(N)
    expected = DenseReference.expect(φ, paulistr(inv(G) * pauliop(str)))
    worst = max(worst, abs(real(expectation(ψ, pauliop(str))) - expected))
  end
  @test worst <= EXACT_ATOL
  @test 0 <= k <= N
  @test any(c != one(CliffordOperator, N) for c in cliffords)    # the Cliffords are not trivial

  # Entropies of the result: SRE is frame-invariant, so it must match the dense state |φ⟩
  # (G_t is a Clifford). This is what Randrot_complexity.jl measures.
  @test sEntropy(ψ; α = 2.0) ≈ dense_sre(φ, N, 2) atol = 1e-9
  @test length(eEntropys!(copy(ψ.mps))) == N - 1
end
