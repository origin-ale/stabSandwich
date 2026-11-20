# T8: the pp → MagicTensors `QubitPauliSum` converter and the `expectation` seam
# (CampsPP/src/MagicTensorsPauliPropagation.jl). Expected values: the bit rule of
# CONVENTIONS.md §3 (helpers.jl `termstring`), the dense reference, and the old-code ⟨O⟩
# saved in data/t8_handoff during T-pre. (The live old-stack comparison on the same circuits
# went with the old code in Phase 6.)

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

  # -- Item 3: converter across the PP ≥ 0.9 key types -----------------------------------
  # Uses `BOUNDARY_N`, `ppstring` and `random_boundary_string` from TestPPKeys.jl.
  # Strings built by hand and passed to PP by symbol (`ppstring`, TestPPKeys.jl), so the
  # expected strings do not go through any bit reading. N covers both sides of each type
  # boundary of the ≥ 0.9 rule (CONVENTIONS.md §3) and the 64-bit limb boundaries above 64
  # qubits (`random_boundary_string` puts Paulis on qubits 32|33, 64|65, …). Runs on any PP
  # version.

  @testset "item 3: converter round trip, N = $N, $T" for N in BOUNDARY_N, T in (Float64, ComplexF64)
    rng = Xoshiro(hash((:item3conv, N, T)))
    rows = Dict{String, T}()
    while length(rows) < min(40, big(4)^N - 1)
      str = N <= 3 ? String(rand(rng, collect("IXYZ"), N)) : random_boundary_string(rng, N)
      all(==('I'), str) && continue
      rows[str] = T <: Complex ? complex(rand(rng, (-1, 1)) * (0.1 + rand(rng)), randn(rng)) :
                                 rand(rng, (-1, 1)) * (0.1 + rand(rng))     # signed, nonzero
    end
    O = pp.PauliSum(T, N)
    for (str, c) in rows
      pp.add!(O, ppstring(str, c))
    end
    @test length(O) == length(rows)
    S = mt.QubitPauliSum(O)
    @test mt.nsites(S) == N
    @test stringdict(S) == Dict(str => ComplexF64(c) for (str, c) in rows)
    @test allunique(objectid.(collect(keys(S.paulis))))
    for (str, c) in Iterators.take(rows, 10)
      t = ppstring(str, c).term
      @test opstring(pauliop_from_ppterm(t, N)) == str
      @test stringdict(mt.QubitPauliSum(ppstring(str, c))) == Dict(str => ComplexF64(c))
    end
  end

  @testset "item 3: empty sums, N = $N" for N in (1, 10, 46, 70)
    S = mt.QubitPauliSum(pp.PauliSum(N))
    @test length(S) == 0
    @test mt.nsites(S) == N
  end

  @testset "item 3: no key aliasing, all 4^3 strings" begin
    strs = [String([a, b, c]) for a in "IXYZ", b in "IXYZ", c in "IXYZ"][:]
    O = pp.PauliSum(Float64, 3)
    for (i, str) in enumerate(strs)
      pp.add!(O, ppstring(str, Float64(i)))
    end
    @test length(O) == 64
    @test stringdict(mt.QubitPauliSum(O)) == Dict(str => ComplexF64(i) for (i, str) in enumerate(strs))
  end

  # -- Seam ------------------------------------------------------------------------------

  @testset "expectation(QubitCAMPS, pp.PauliSum) vs dense, N = $N" for N in (3, 4, 5)
    rng = Xoshiro(hash((:t8seam, N)))
    for _ in 1:4
      ones = [j for j in 1:N if rand(rng, Bool)]
      gates, phases = rotation_circuit(rng, 8, N)
      # Mix in Clifford angles, so the frame is a non-trivial Clifford (the first is pinned
      # to +π/4 so the frame is surely non-Pauli).
      cl = rand(rng, eachindex(phases), 3)
      phases[cl] .= rand(rng, (π/4, -π/4, π/2, -π/2, 3π/4, -3π/4, 0.0, π), 3)
      phases[cl[1]] = π/4
      dense, ψ = both_states(mt.QubitCAMPS, N, ones, gates, phases)
      @test mt.get_clifford_copy(ψ).clifford != one(CliffordOperator, N)
      for T in (Float64, ComplexF64)
        O = random_ppsum(rng, N, min(40, 4^N - 1); T)
        ev = mt.expectation(ψ, O)
        @test ev ≈ dense_expect(dense, O) atol = EXACT_ATOL
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
    # Rebuild with the reference's own stopping rule (generate.jl: CAMPS until the bond
    # dimension reaches χ = 64, which first happens at gate s). A plain `evolve` to gate s has
    # no bond-dimension stop, so a broken frame or disentangler made the bond dimension blow up
    # and hung the suite (T10 mutants E7, E8, D5); this rebuild stops as soon as it reaches 64.
    ψ, k, s_rebuilt = DisentangleCAMPS.evolve_bonddim(
      computationalcamps(mt.QubitCAMPS, N, excited), 64, gates, phases)
    @test s_rebuilt == s
    # Frame-dependent parity with the old code (state.txt): k and the bond dimension.
    @test k == 37
    @test DisentangleCAMPS.bonddim(ψ) == 112
    # The sums only mean something on the reference state; skip them on a wrong rebuild
    # (its failures are already reported above).
    files = s_rebuilt == s ? (SLOW ? ("sum_0201.txt", "sum_0538.txt") : ("sum_0201.txt",)) : ()
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
      # Post-port item 4 (CONVENTIONS.md §5): the per-term default and the kept MPO route
      # both reproduce the old value and agree to 1e-12.
      evp = mt.expectation(ψ, O; method = :perterm)
      evm = mt.expectation(ψ, O; method = :mpo)
      @test abs(evp - old) <= EXACT_ATOL
      @test abs(evm - old) <= EXACT_ATOL
      @test abs(evp - evm) <= 1e-12
    end
  end
end
