# Initialize.jl against CONVENTIONS.md §1 and §5.

"⟨Z_j⟩ of a CAMPS (either stack)."
zexp(ψ, j) = real(state_expectation(ψ, single_z(nqubits_state(ψ), j)))

@testset "Initialize" begin
  @testset "computationalcamps, $(stackname(T))" for T in STATE_TYPES
    for (N, ones) in [(3, Int[]), (4, [1, 4]), (5, [2, 3, 5]), (6, collect(1:6))]
      ψ = computationalcamps(T, N, ones)
      @test ψ isa T
      @test nqubits_state(ψ) == N
      @test all(zexp(ψ, j) ≈ (j in ones ? -1.0 : 1.0) for j in 1:N)   # |1⟩ ↔ ⟨Z⟩ = −1
      @test real(state_expectation(ψ, single_x(N, 1))) ≈ 0 atol = EXACT_ATOL
      @test DisentangleCAMPS.bonddim(ψ) == 1                          # frame-dependent
    end
    # The untyped form still builds the old type while the port is staged (plan, Phase 2).
    @test computationalcamps(3, [1]) isa cmps.CAMPS
  end

  @testset "domainwallstate, $(stackname(T))" for T in STATE_TYPES
    # The state matches the returned list, and the draws don't depend on the state type.
    N = 10
    ψ, ones = domainwallstate(T, Xoshiro(3), N, 0.6)
    @test ψ isa T
    @test all(zexp(ψ, j) ≈ (j in ones ? -1.0 : 1.0) for j in 1:N)
    @test ones == domainwallstate(Xoshiro(3), N, 0.6)[2]
  end

  @testset "domainwallstate" begin
    N = 10
    ψ, ones = domainwallstate(Xoshiro(3), N, 0.6)
    @test issorted(ones) && allunique(ones) && all(1 .<= ones .<= N)
    @test all(zexp(ψ, j) ≈ (j in ones ? -1.0 : 1.0) for j in 1:N)     # state ↔ returned list
    @test domainwallstate(Xoshiro(3), N, 0.6)[2] == ones               # reproducible

    # μ = 0: every qubit is |1⟩ with probability 1/2, on both halves. Fixed seeds, 5σ.
    n = 400
    counts_l = 0; counts_r = 0
    for seed in 1:n
      _, o = domainwallstate(Xoshiro(seed), N, 0.0)
      counts_l += count(<=(N ÷ 2), o); counts_r += count(>(N ÷ 2), o)
    end
    σ = sqrt(n * (N ÷ 2) * 0.25)
    @test abs(counts_l - n * (N ÷ 2) / 2) <= 5σ
    @test abs(counts_r - n * (N ÷ 2) / 2) <= 5σ

    # μ = 0.6 (the production value): p = e^0.6 / (e^0.6 + e^−0.6) = 0.7685. Each left qubit
    # is |1⟩ with probability 1 − p, each right qubit with probability p, independently.
    # Per-qubit frequencies over 400 fixed seeds, 5σ (σ ≈ 0.021); a different formula such
    # as e^μ/(1 + e^μ) = 0.646 is more than 5σ away.
    p = exp(0.6) / (exp(0.6) + exp(-0.6))
    freq = zeros(N)
    for seed in 1:n
      _, o = domainwallstate(Xoshiro(seed), N, 0.6)
      freq[o] .+= 1 / n
    end
    σq = sqrt(p * (1 - p) / n)
    @test all(abs.(freq[1:N÷2] .- (1 - p)) .<= 5σq)
    @test all(abs.(freq[N÷2+1:N] .- p) .<= 5σq)
    @test abs(p - 0.646) > 5σq

    # Orientation (CONVENTIONS.md §1): with p = e^μ / (e^μ + e^−μ), a left qubit is |1⟩ with
    # probability 1 − p and a right qubit with probability p, so for large μ the
    # excitations sit on the right half.
    _, o = domainwallstate(Xoshiro(5), N, 20.0)
    @test o == collect(N÷2+1:N)
    _, o = domainwallstate(Xoshiro(5), N, -20.0)
    @test o == collect(1:N÷2)
  end

  @testset "transferredmagnetization: operator, N=$N" for (N, init) in
      [(4, [3, 4]), (5, [1, 4]), (6, [1, 4]), (6, [4, 5, 6]), (7, Int[]), (46, collect(24:46))]
    # O = −½ Σ_{j∈L} Z_j + ½ Σ_{j∈R} Z_j + (N/2 − n₁ᴸ(0) − n₀ᴿ(0))·1, L = 1…⌊N/2⌋.
    O = transferredmagnetization(N, init)
    L = 1:fld(N, 2)
    n1L = count(in(L), init)
    n0R = count(j -> !(j in L) && !(j in init), 1:N)
    T = keytype(O.terms)
    zterm(j) = T(3) << (2 * (j - 1))               # Z on qubit j, by the §3 bit rule

    @test O isa pp.PauliSum
    @test length(O.terms) == N + 1                 # N Z terms and the identity, always
    @test all(O.terms[zterm(j)] == (j in L ? -0.5 : 0.5) for j in 1:N)
    @test O.terms[zero(T)] == N / 2 - n1L - n0R

    # ⟨O⟩ = 0 on the initial state.
    @test pp.overlapwithcomputational(O, init) ≈ 0 atol = EXACT_ATOL
    if N <= 7
      @test dense_expect(DR.basis_state(N, init), O) ≈ 0 atol = EXACT_ATOL
    end
  end

  @testset "transferredmagnetization: values on basis states" begin
    N = 6
    val(init, final) = dense_expect(DR.basis_state(N, final), transferredmagnetization(N, init))
    @test val([4, 5], [3, 5]) ≈ 2        # one excitation crossed right → left
    @test val([4, 5], [1, 2]) ≈ 4        # two crossed
    @test val([2, 4, 5], [4, 5, 6]) ≈ -2 # one crossed left → right
    @test val([4, 5], [3, 4, 5]) ≈ 1     # a 1 created on the left (ΣZ not conserved)
    @test val([4, 5], [4]) ≈ 1           # a 1 removed on the right
    @test val([4, 5], [5, 4]) ≈ 0        # nothing moved

    # The same numbers through the code's own evaluation paths.
    O = transferredmagnetization(N, [4, 5])
    @test pp.overlapwithcomputational(O, [3, 5]) ≈ 2
    for T in STATE_TYPES
      @test real(state_expectation(computationalcamps(T, N, [3, 5]), O)) ≈ 2 atol = EXACT_ATOL
      @test real(state_expectation(computationalcamps(T, N, [1, 2]), O)) ≈ 4 atol = EXACT_ATOL
    end
  end
end
