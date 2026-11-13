# T8: the pp → MagicTensors `QubitPauliSum` converter and the `expectation` seam
# (CampsPP/src/MagicTensorsPauliPropagation.jl). Expected values: the bit rule of
# CONVENTIONS.md §3 (helpers.jl `termstring`), the dense reference, the old CliffordMPS code
# on the same physical state, and the old-code ⟨O⟩ saved in data/t8_handoff during T-pre.

using ITensors, ITensorMPS
import MagicTensors as mt

const SLOW = get(ENV, "STABSANDWICH_SLOW", "") != ""

# -- Helpers -----------------------------------------------------------------------------

"Character of qubit j of a QuantumClifford Pauli string (I, X, Y, Z), sign ignored."
paulichar(P::PauliOperator, j) = "IXZY"[1 + xbit(P)[j] + 2 * zbit(P)[j]]
opstring(P::PauliOperator) = String([paulichar(P, j) for j in 1:nqubits(P)])

"A converted sum as Dict(string => coefficient); the stored strings carry no sign."
function stringdict(S::mt.QubitPauliSum)
  d = Dict{String, ComplexF64}()
  for (Q, c) in S
    @assert iszero(Q.pauli.phase[])
    d[opstring(Q.pauli)] = c
  end
  return d
end

"The expected Dict(string => coefficient) of a pp sum, by the §3 bit rule."
spec_stringdict(O::pp.PauliSum) =
  Dict(termstring(t, O.nqubits) => ComplexF64(c) for (t, c) in O.terms)

"pp sum from (string, coefficient) rows, character j = qubit j (as generate.jl)."
function ppsum_from_rows(N, rows; T = Float64)
  S = pp.PauliSum(T, N)                    # pp.PauliSum(N) would force Float64 coefficients
  for (str, c) in rows
    qs = [j for j in 1:N if str[j] != 'I']
    pp.add!(S, pp.PauliString(N, [Symbol(str[j]) for j in qs], qs, T(c)))
  end
  return S
end

"Random pp sum of about n distinct terms on N qubits, with terms on qubit N and qubits > 32."
function random_ppsum(rng, N, n; T = Float64)
  n <= big(4)^N - 1 || throw(ArgumentError("only 4^N − 1 = $(big(4)^N - 1) non-identity strings on N = $N qubits, asked for $n"))
  rows = Dict{String, T}()
  while length(rows) < n
    w = rand(rng, 1:min(N, 6))
    qs = unique(vcat(rand(rng, 1:N, w), rand(rng) < 0.3 ? [N] : Int[]))
    chars = fill('I', N)
    for q in qs
      chars[q] = rand(rng, "XYZ")
    end
    rows[String(chars)] = T <: Complex ? randn(rng, T) : randn(rng)
  end
  return ppsum_from_rows(N, collect(rows); T)
end

"State after the circuit, on the dense reference and on a CAMPS of type T."
function both_states(T, N, ones, gates, phases)
  dense = DR.basis_state(N, ones)
  for (g, ϕ) in zip(gates, phases)
    DR.rotate!(dense, gatestring(g, N), ϕ)
  end
  ψ, _ = evolve(computationalcamps(T, N, ones), length(gates), gates, phases)
  return dense, ψ
end

# -- Converter ---------------------------------------------------------------------------

@testset "T8: pp → QubitPauliSum converter and expectation seam" begin
  @testset "hand cases (CONVENTIONS.md §3)" begin
    S = mt.QubitPauliSum(pp.PauliString(3, [:X, :Y, :Z], [1, 2, 3], -3.0))
    @test stringdict(S) == Dict("XYZ" => -3.0 + 0im)
    S = mt.QubitPauliSum(pp.PauliString(46, [:Y, :Z], [33, 46], 0.25))
    @test stringdict(S) == Dict("I"^32 * "Y" * "I"^12 * "Z" => 0.25 + 0im)
    S = mt.QubitPauliSum(pp.PauliSum(2, Dict(UInt8(9) => 1.0, UInt8(3) => 0.5)))  # XY, ZI
    @test stringdict(S) == Dict("XY" => 1.0 + 0im, "ZI" => 0.5 + 0im)
    S = mt.QubitPauliSum(pp.PauliString(2, [:I, :I], [1, 2], 0.5))
    @test stringdict(S) == Dict("II" => 0.5 + 0im)
    # Coefficients become ComplexF64, whatever the pp coefficient type.
    @test valtype(mt.QubitPauliSum(pp.PauliString(4, [:X], [2], 1.5)).paulis) == ComplexF64
    @test valtype(mt.QubitPauliSum(pp.PauliString(4, [:X], [2], 1.5 - 2im)).paulis) == ComplexF64
    @test stringdict(mt.QubitPauliSum(pp.PauliString(4, [:X], [2], 1.5 - 2im))) == Dict("IXII" => 1.5 - 2im)
  end

  @testset "random sums against the bit rule, N = $N, $T" for N in (1, 4, 33, 46, 70),
      T in (Float64, ComplexF64)
    rng = Xoshiro(hash((:t8conv, N, T)))
    n = N == 1 ? 3 : Int(min(300, (big(4)^N - 1) ÷ 2))   # distinct strings must exist (4^N − 1 of them)
    O = random_ppsum(rng, N, n; T)
    S = mt.QubitPauliSum(O)
    @test mt.nsites(S) == N
    @test length(S) == length(O)
    @test stringdict(S) == spec_stringdict(O)                  # strings, signs and values
    # A fresh key vector per term: no two terms share storage (a shared buffer would make
    # every key alias the last one written).
    @test allunique(objectid.(collect(keys(S.paulis))))
    # Same, through the single-string converter and `pauliop_from_ppterm`.
    for (t, c) in Iterators.take(O.terms, 20)
      @test stringdict(mt.QubitPauliSum(pp.PauliString(N, t, c))) == Dict(termstring(t, N) => ComplexF64(c))
      @test opstring(pauliop_from_ppterm(t, N)) == termstring(t, N)
    end
    # Sensitivity: a wrong qubit order would be noticed on these sums.
    N > 1 && @test any(reverse(k) != k for k in keys(spec_stringdict(O)))
  end

  @testset "converting does not change the pp sum" begin
    O = random_ppsum(Xoshiro(3), 46, 50)
    before = copy(O.terms)
    mt.QubitPauliSum(O)
    @test O.terms == before
  end

  # -- Seam ------------------------------------------------------------------------------

  @testset "expectation(QubitCAMPS, pp.PauliSum) vs dense and old code, N = $N" for N in (3, 4, 5)
    rng = Xoshiro(hash((:t8seam, N)))
    for trial in 1:4
      ones = [j for j in 1:N if rand(rng, Bool)]
      gates, phases = rotation_circuit(rng, 8, N)
      # Mix in Clifford angles, so the frame is a non-trivial Clifford.
      cl = rand(rng, eachindex(phases), 3)
      phases[cl] .= rand(rng, (π/4, -π/4, π/2, 3π/4), 3)
      dense, ψ = both_states(mt.QubitCAMPS, N, ones, gates, phases)
      @test mt.get_clifford_copy(ψ).clifford != one(CliffordOperator, N)
      for T in (Float64, ComplexF64)
        O = random_ppsum(rng, N, min(40, 4^N - 1); T)
        ev = mt.expectation(ψ, O)
        @test ev ≈ dense_expect(dense, O) atol = EXACT_ATOL
        # The old stack on the same circuit: same physical value. (Old code only on angles
        # it handles correctly: no −π/4, +3π/4 or ±π/2 there, hard spot 6.)
        if all(ϕ -> !(ϕ in (-π/4, 3π/4, π/2)), phases)
          ψold, _ = evolve(computationalcamps(cmps.CAMPS, N, ones), length(gates), gates, phases)
          @test ev ≈ cmps.expectation(ψold, O) atol = EXACT_ATOL
        end
      end
      # Signed single strings (PauliOperator method): the sign counts.
      for _ in 1:10
        P = random_pauli(rng, N; nophase = true)
        s = rand(rng, (1, -1))
        Ps = s == 1 ? P : -P
        @test mt.expectation(ψ, Ps) ≈ s * DR.expect(dense, opstring(P)) atol = EXACT_ATOL
      end
    end
  end

  @testset "verbose = true gives the same value" begin
    N = 4
    dense, ψ = both_states(mt.QubitCAMPS, N, [2], rotation_circuit(Xoshiro(5), 6, N)...)
    O = random_ppsum(Xoshiro(6), N, 30)
    ev = mt.expectation(ψ, O)
    ev_verbose = redirect_stdout(devnull) do
      mt.expectation(ψ, O; verbose = true)
    end
    @test ev_verbose == ev
  end

  # -- Real N = 46 hand-off sums (old-code reference) --------------------------------------

  @testset "N = 46 hand-off sums reproduce the old ⟨O⟩" begin
    dir = joinpath(@__DIR__, "data", "t8_handoff")
    fields = Dict(first(split(l, '\t')) => last(split(l, '\t'; limit = 2))
                  for l in eachline(joinpath(dir, "state.txt")) if !startswith(l, '#'))
    N = parse(Int, fields["N"])
    s = parse(Int, fields["s"])
    excited = parse.(Int, split(fields["excited"]))
    doped = parse.(Int, split(fields["doped_gates"]))
    gates, phases = xxz_circuit(π/4, π/4, 23, N)
    phases[doped] .= 3π/16
    ψ, k = evolve(computationalcamps(mt.QubitCAMPS, N, excited), s, gates, phases)
    # Frame-dependent parity with the old code (state.txt): k and the bond dimension.
    @test k == 37
    @test DisentangleCAMPS.bonddim(ψ) == 112
    files = SLOW ? ("sum_0201.txt", "sum_0538.txt") : ("sum_0201.txt",)
    @testset "$f" for f in files
      lines = collect(eachline(joinpath(dir, f)))
      header = Dict(split(l[3:end], '\t')[1] => split(l[3:end], '\t')[2]
                    for l in lines if startswith(l, "# ") && occursin('\t', l))
      rows = [(first(split(l, '\t')), parse(Float64, last(split(l, '\t'))))
              for l in lines if !startswith(l, '#')]
      O = ppsum_from_rows(N, rows)
      @test length(O) == parse(Int, header["nterms"])
      old = parse(ComplexF64, replace(header["old_expectation_mpo"], " " => ""))
      ev = mt.expectation(ψ, O)
      @test abs(ev - old) <= EXACT_ATOL
      @test abs(imag(ev)) <= 1e-12
    end
  end
end
