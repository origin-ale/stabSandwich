# Post-port item 4: the pp seam
#   expectation(ψ::QubitCAMPS, O::pp.PauliSum; method = :perterm, normalized = false, verbose = false)
# (test/CONVENTIONS.md §5, "The pp seam"; approved 2026-11-19). Written from the spec before the
# per-term code (T8-style). Expected values: hand derivation (comments), the dense reference
# (N ≤ 8) and the stabilizer reference (N = 46). The old-code values of the T8 hand-off sums
# are checked in TestMagicTensorsPP.jl. Nothing here is read off the implementation.
#
# Gating. `SEAM_ITEM4` (helpers.jl) is true once the seam accepts the `method` and `normalized` keywords.
# Until then:
#   - the groups marked [item 4] record one `@test_broken false` each instead of running.
#     That covers the new keywords, and also two spec points the pre-item-4 seam does not meet:
#     an empty sum (the MPO route throws) and ψ left bitwise unchanged (upstream's
#     `get_clifford(ψ; daggered = true)` inverts a stored C in place and flips `daggered`);
#   - the rest (values with the default keywords, hand pins, the reference self-checks) runs on
#     any seam.
# When the code lands, the gate opens by itself and the broken markers disappear.
#
# States. The CAMPS states are built directly from their two parts, so the frame C and the MPS
# are known without asking the state: a random MPS (random complex tensors, bond dimensions
# 1–8, not normalised) in a chosen gauge, and a random Clifford C stored as C or as C†
# (CONVENTIONS.md §4). The reference value is Σ c ⟨mps|C†PC|mps⟩, with C†PC from
# QuantumClifford's `C⁻¹ * Stabilizer` (pinned as conjugation in TestGoldenSet.jl) and
# ⟨mps|Q|mps⟩ from the dense reference on the contracted MPS.
#
# Sums. A per-term evaluation walks the conjugated strings Q = C†PC, so the structures the plan
# asks for (shared prefixes, strings differing only at their last site, last non-identity site
# 1 or N, the identity) are built on Q and mapped back to P = C Q C† (`sum_from_Q`).

using ITensors, ITensorMPS
using LinearAlgebra: norm
import MagicTensors as mt

isdefined(@__MODULE__, :StabilizerReference) ||
  include(joinpath(@__DIR__, "DisentangleCAMPS", "StabilizerReference.jl"))
isdefined(@__MODULE__, :SR) || @eval const SR = StabilizerReference   # also set by TestCliffordLargeN.jl

const SEAM_SLOW = get(ENV, "STABSANDWICH_SLOW", "") != ""


# -- Strings and signs ---------------------------------------------------------------------

"Character of qubit j of a QuantumClifford Pauli string (I, X, Y, Z), sign ignored."
seam_char(P::PauliOperator, j) = "IXZY"[1 + xbit(P)[j] + 2 * zbit(P)[j]]
seam_opstring(P::PauliOperator) = String([seam_char(P, j) for j in 1:nqubits(P)])

"String over IXYZ (character j = qubit j) → Hermitian PauliOperator with sign +."
seam_strop(str) = PauliOperator(0x0, [c in "XY" for c in str], [c in "ZY" for c in str])

"Sign of a Hermitian PauliOperator (§3: phase byte 0 → +1, 2 → −1)."
function seam_sign(P::PauliOperator)
  P.phase[] == 0x0 && return 1
  P.phase[] == 0x2 && return -1
  error("non-Hermitian string: bug in the test")
end

"A P A† for a CliffordOperator A (QC `A * Stabilizer` is conjugation, CONVENTIONS.md §4)."
seam_conj(A::CliffordOperator, P::PauliOperator) = (A * Stabilizer([P]))[1]

"pp sum from (string, coefficient) rows; the all-I string becomes the identity term."
function seam_ppsum(N, rows; T = Float64)
  S = pp.PauliSum(T, N)
  for (str, c) in rows
    qs = [j for j in 1:N if str[j] != 'I']
    syms = [Symbol(str[j]) for j in qs]
    isempty(qs) && (qs = [1]; syms = [:I])
    pp.add!(S, pp.PauliString(N, syms, qs, T(c)))
  end
  return S
end

# -- States ----------------------------------------------------------------------------------

"""
Random MPS on N "Qubit" sites: link dimensions drawn from 1:χmax, i.i.d. complex Gaussian
tensors, no orthogonality centre, rescaled so that ⟨mps|mps⟩ = `norm2` (not 1 by default, so
that normalisation is visible).
"""
function random_raw_mps(rng, N, χmax; norm2 = 0.5 + 1.5 * rand(rng))
  s = siteinds("Qubit", N)
  l = [Index(rand(rng, 1:χmax), "Link,l=$b") for b in 1:N-1]
  tensors = map(1:N) do j
    is = Index[]
    j > 1 && push!(is, l[j-1])
    push!(is, s[j])
    j < N && push!(is, l[j])
    ITensor(randn(rng, ComplexF64, dim.(is)...), is...)
  end
  m = MPS(tensors)
  m[1] *= sqrt(norm2) / norm(m)
  return m
end

"The MPS in gauge `g`: :raw (no centre), or (:centre, j) after orthogonalize! to site j."
function in_gauge(m, g)
  g === :raw && return deepcopy(m)
  return orthogonalize!(deepcopy(m), g[2])
end

"""
Dense vector of an MPS in the layout of the dense reference (§1: the amplitude of |b₁…b_N⟩ at
1 + Σⱼ bⱼ 2^(j−1)). A "Qubit" site has value 1 for |0⟩, 2 for |1⟩; Julia arrays are column-major,
so listing the sites as s₁, …, s_N makes qubit 1 the fastest index. Pinned below.
"""
function mps_dense(m)
  s = siteinds(m)
  return vec(Array(foldl(*, [m[j] for j in eachindex(m)]), s...))
end

"""
CAMPS with state frame C (|ψ⟩ = C|mps⟩, §4), stored as C (`daggered = false`) or as C†
(`daggered = true`). Both parts are copied.
"""
function camps_from(m, C::CliffordOperator; daggered::Bool)
  stored = daggered ? inv(C) : copy(C)
  return mt.QubitCAMPS(deepcopy(m), mt.QubitCliffordUnitary(stored), daggered)
end

# -- Reference values ------------------------------------------------------------------------

"""
Spec value (§5): Σ c ⟨mps|C†PC|mps⟩ over the terms (P, c) of O, with v the dense MPS and
Cinv = C⁻¹ (so C†PC = Cinv P Cinv†). The identity term gives c·⟨v|v⟩, an empty sum 0.
"""
function spec_seam(v, Cinv::CliffordOperator, O::pp.PauliSum)
  N = O.nqubits
  ev = zero(ComplexF64)
  for (t, c) in O.terms
    Q = seam_conj(Cinv, seam_strop(termstring(t, N)))
    ev += ComplexF64(c) * seam_sign(Q) * DR.expect(v, seam_opstring(Q))
  end
  return ev
end

# -- Sums designed on the conjugated strings --------------------------------------------------

"Random string over IXYZ of length n."
randstr(rng, n) = String(rand(rng, collect("IXYZ"), n))
"Random non-identity letter."
randnonI(rng) = rand(rng, ('X', 'Y', 'Z'))

"""
Conjugated strings Q (unsigned, character j = qubit j) with the structures of the plan:
the identity; last non-identity site 1; last non-identity site N; families sharing a long
prefix; pairs equal up to their last non-identity site m and different there (that is,
differing only at the last site); and plain random strings. About `n` distinct strings.
"""
function designed_Qs(rng, N, n)
  Qs = Set{String}()
  push!(Qs, "I"^N)                                             # identity
  for a in "XYZ"                                               # last site 1
    push!(Qs, a * "I"^(N - 1))
  end
  for _ in 1:3                                                 # last site N
    push!(Qs, randstr(rng, N - 1) * randnonI(rng))
  end
  while length(Qs) < n
    r = rand(rng)
    if r < 0.4                                                 # shared long prefix, varied tails
      L = rand(rng, max(1, N - 3):N-1)
      prefix = randstr(rng, L)
      for _ in 1:rand(rng, 2:5)
        push!(Qs, prefix * randstr(rng, N - L - 1) * randnonI(rng))
      end
    elseif r < 0.7                                             # differ only at the last site m
      m = rand(rng, 1:N)
      head = randstr(rng, m - 1)
      for a in "XYZ"
        push!(Qs, head * a * "I"^(N - m))
      end
    else
      push!(Qs, randstr(rng, N))
    end
  end
  return collect(Qs)
end

"Random coefficient of type T, bounded away from 0, with both signs."
randcoeff(rng, ::Type{Float64}) = rand(rng, (-1, 1)) * (0.1 + rand(rng))
randcoeff(rng, ::Type{ComplexF64}) = complex(randcoeff(rng, Float64), randcoeff(rng, Float64))

"""
pp sum whose conjugated strings through the frame C are `Qs`: each Q becomes the unsigned
string of P = C Q C† with a random coefficient. Conjugation is a bijection on strings, so the
terms stay distinct.
"""
function sum_from_Q(rng, C::CliffordOperator, Qs; T = Float64)
  rows = [(seam_opstring(seam_conj(C, seam_strop(Q))), randcoeff(rng, T)) for Q in Qs]
  @assert allunique(first.(rows))
  return seam_ppsum(nqubits(C), rows; T)
end

# -- Large-N states ---------------------------------------------------------------------------

# The MPS is a stabilizer state made with ITensors gates (H, S, CNOT at any distance), so
# its bond dimension is above 1; its generators follow the same gates through
# QuantumClifford's symbolic gates. Assumption pins: ITensors' "S" is diag(1, i) and "CNOT"
# has its control on the first site; the N = 4 self-check compares the two routes on all
# 256 strings through the dense reference.
"Stabilizer MPS and its generators (a Stabilizer) after `ngates` random H, S, CNOT gates."
function stabilizer_mps(rng, N, ngates; maxdist = 6)
  s = siteinds("Qubit", N)
  m = MPS(ComplexF64, s, "0")
  S = Stabilizer(SR.basis_generators(N, Int[]))
  for _ in 1:ngates
    r = rand(rng)
    if r < 0.3
      q = rand(rng, 1:N)
      m = apply(op("H", s[q]), m)
      apply!(S, sHadamard(q))
    elseif r < 0.5
      q = rand(rng, 1:N)
      m = apply(op("S", s[q]), m)
      apply!(S, sPhase(q))
    else
      c = rand(rng, 1:N)
      t = clamp(c + rand(rng, (-1, 1)) * rand(rng, 1:maxdist), 1, N)
      t == c && continue
      m = apply(op("CNOT", s[c], s[t]), m; cutoff = 1e-13)
      apply!(S, sCNOT(c, t))
    end
  end
  return m, S
end

# -- Snapshots for "ψ is unchanged" -----------------------------------------------------------

"Deep copy of everything the spec says must stay bitwise equal: tensors with indices, storage."
snapshot(ψ::mt.QubitCAMPS) =
  (tensors = [deepcopy(ψ.mps[j]) for j in eachindex(ψ.mps)],
   tableau = deepcopy(ψ.clifford.clifford), daggered = ψ.daggered)

"Is ψ bitwise equal to the snapshot? (Indices as sets and in order; data in the stored order.)"
function unchanged(ψ::mt.QubitCAMPS, snap)
  length(ψ.mps) == length(snap.tensors) || return false
  for (A, B) in zip(ψ.mps, snap.tensors)
    inds(A) == inds(B) || return false
    Array(A, inds(B)...) == Array(B, inds(B)...) || return false
  end
  return ψ.daggered == snap.daggered && ψ.clifford.clifford == snap.tableau
end

"One broken marker for an [item 4] group while the seam lacks the item 4 code."
item4_pending() = @test_broken false

# ============================================================================================

@testset "item 4: pp seam (CONVENTIONS.md §5)" begin
  @info "item 4 seam tests: SEAM_ITEM4 = $SEAM_ITEM4" *
        (SEAM_ITEM4 ? "" : " (groups marked [item 4] are recorded as broken)")

  # -- Self-checks of the test helpers ----------------------------------------------------

  @testset "helpers: MPS → dense layout and the conjugation route" begin
    # |0 1 1 0 0⟩ sits at 1 + 2 + 4 = 7 (§1).
    s = siteinds("Qubit", 5)
    v = mps_dense(MPS(s, ["0", "1", "1", "0", "0"]))
    @test v ≈ DR.basis_state(5, [2, 3])
    # Random MPS: every ⟨P⟩ from the dense route equals the MPO-free ITensor contraction
    # ⟨m|P|m⟩ with ITensors' own Pauli ops (a cross-check of the layout, N = 3).
    m = random_raw_mps(Xoshiro(1), 3, 3)
    v = mps_dense(m)
    sm = siteinds(m)
    for str in DR.all_paulis(3)
      opsym = Dict('I' => "Id", 'X' => "X", 'Y' => "Y", 'Z' => "Z")
      Pm = apply([op(opsym[str[j]], sm[j]) for j in 1:3], m)
      @test DR.expect(v, str) ≈ real(inner(m, Pm)) atol = EXACT_ATOL
    end
    # Negative check: the reversed layout gives other values on this state.
    vrev = vec(permutedims(reshape(v, 2, 2, 2), (3, 2, 1)))
    @test any(abs(DR.expect(vrev, str) - DR.expect(v, str)) > 1e-3 for str in DR.all_paulis(3))
    # `designed_Qs` builds what it promises.
    Qs = designed_Qs(Xoshiro(2), 6, 60)
    @test "IIIIII" in Qs
    @test any(q -> q[1] != 'I' && all(==('I'), q[2:end]), Qs)
    @test any(q -> q[end] != 'I', Qs)
    lastsite(q) = something(findlast(!=('I'), q), 0)
    @test any(q1 != q2 && lastsite(q1) == lastsite(q2) > 1 &&
              q1[1:lastsite(q1)-1] == q2[1:lastsite(q2)-1] for q1 in Qs, q2 in Qs)
    # `sum_from_Q`: conjugating each P back through C† gives ±Q.
    C = random_clifford(Xoshiro(3), 6)
    O = sum_from_Q(Xoshiro(4), C, Qs)
    @test length(O) == length(Qs)
    @test Set(seam_opstring(seam_conj(inv(C), seam_strop(termstring(t, 6)))) for t in keys(O.terms)) == Set(Qs)
  end

  # -- Hand-derived pins (default keywords) ------------------------------------------------

  @testset "hand pins" begin
    # N = 2, identity frame, |mps⟩ = |01⟩: ⟨Z₁⟩ = +1, ⟨Z₂⟩ = −1, ⟨X₁⟩ = 0, ⟨II⟩ = 1.
    s = siteinds("Qubit", 2)
    ψ = camps_from(MPS(s, ["0", "1"]), one(CliffordOperator, 2); daggered = false)
    O = seam_ppsum(2, [("ZI", 1.0), ("IZ", 10.0), ("XI", 100.0), ("II", 1000.0)])
    @test mt.expectation(ψ, O) ≈ 1 - 10 + 0 + 1000 atol = EXACT_ATOL
    # C = H₁: |ψ⟩ = H|0⟩ ⊗ |0⟩ = |+0⟩, so ⟨X₁⟩ = 1, ⟨Z₁⟩ = 0, ⟨Z₂⟩ = 1.
    H1 = CliffordOperator(T"Z_ _X X_ _Z")          # X₁ ↦ Z₁, X₂ ↦ X₂, Z₁ ↦ X₁, Z₂ ↦ Z₂
    for d in (false, true)
      ψ = camps_from(MPS(s, "0"), H1; daggered = d)
      @test mt.expectation(ψ, seam_ppsum(2, [("XI", 1.0)])) ≈ 1 atol = EXACT_ATOL
      @test mt.expectation(ψ, seam_ppsum(2, [("ZI", 1.0)])) ≈ 0 atol = EXACT_ATOL
      @test mt.expectation(ψ, seam_ppsum(2, [("IZ", 2.0), ("XI", -0.5)])) ≈ 1.5 atol = EXACT_ATOL
    end
    # C = Z₁·H₁ (H acts first, §4): |ψ⟩ = Z|+⟩ = |−⟩, so ⟨X₁⟩ = −1 (the frame sign).
    Z1 = CliffordOperator(T"-X_ _X Z_ _Z")
    for d in (false, true)
      ψ = camps_from(MPS(s, "0"), Z1 * H1; daggered = d)
      @test mt.expectation(ψ, seam_ppsum(2, [("XI", 1.0)])) ≈ -1 atol = EXACT_ATOL
    end
    # C = S₁·H₁ with S = diag(1, i) (QC tPhase: X ↦ Y, §4): |ψ⟩ = S|+⟩ = |+i⟩, ⟨Y₁⟩ = +1.
    S1 = CliffordOperator(T"Y_ _X Z_ _Z")
    for d in (false, true)
      ψ = camps_from(MPS(s, "0"), S1 * H1; daggered = d)
      @test mt.expectation(ψ, seam_ppsum(2, [("YI", 1.0)])) ≈ 1 atol = EXACT_ATOL
    end
    # Complex coefficient on |+0⟩: (2 − 3i)·⟨X₁⟩ = 2 − 3i; the imaginary part counts.
    ψ = camps_from(MPS(s, "0"), H1; daggered = false)
    ev = mt.expectation(ψ, seam_ppsum(2, [("XI", 2.0 - 3.0im)]; T = ComplexF64))
    @test ev isa ComplexF64
    @test abs(ev - (2 - 3im)) <= EXACT_ATOL
    # Not normalised (§5): |mps⟩ = 2|00⟩ gives ⟨II⟩ = 4 and ⟨Z₁⟩ = 4.
    m = MPS(s, "0")
    m[1] *= 2
    ψ = camps_from(m, one(CliffordOperator, 2); daggered = false)
    @test mt.expectation(ψ, seam_ppsum(2, [("II", 1.0)])) ≈ 4 atol = EXACT_ATOL
    @test mt.expectation(ψ, seam_ppsum(2, [("ZI", 1.0)])) ≈ 4 atol = EXACT_ATOL
  end

  # -- Randomized-exact against the dense reference (default keywords) -------------------

  ntrials = SEAM_SLOW ? 6 : 2
  @testset "random CAMPS vs dense, N = $N" for N in 2:8
    rng = Xoshiro(hash((:item4dense, N)))
    for trial in 1:ntrials, d in (false, true)
      χmax = rand(rng, 1:8)
      m0 = random_raw_mps(rng, N, χmax)
      v = mps_dense(m0)
      C = random_clifford(rng, N)
      Cinv = inv(C)
      for g in (:raw, (:centre, 1), (:centre, N), (:centre, rand(rng, 1:N)))
        ψ = camps_from(in_gauge(m0, g), C; daggered = d)
        for T in (Float64, ComplexF64)
          O = sum_from_Q(rng, C, designed_Qs(rng, N, min(40, 4^N)); T)
          ev = mt.expectation(ψ, O)
          @test ev isa ComplexF64
          @test abs(ev - spec_seam(v, Cinv, O)) <= EXACT_ATOL
        end
      end
    end
  end

  # -- Large N: Clifford MPS behind a random frame, stabilizer reference ------------------

  @testset "assumption pins: ITensors gate matrices" begin
    s = siteinds("Qubit", 2)
    @test Array(op("S", s[1]), s[1]', s[1]) ≈ [1 0; 0 im]
    @test Array(op("H", s[1]), s[1]', s[1]) ≈ [1 1; 1 -1] / sqrt(2)
    # CNOT(control s₁, target s₂): |10⟩ ↦ |11⟩, |00⟩ ↦ |00⟩.
    m = apply(op("CNOT", s[1], s[2]), MPS(ComplexF64, s, ["1", "0"]))
    @test mps_dense(m) ≈ DR.basis_state(2, [1, 2])
    m = apply(op("CNOT", s[1], s[2]), MPS(ComplexF64, s, ["0", "1"]))
    @test mps_dense(m) ≈ DR.basis_state(2, [2])
  end

  @testset "self-check: stabilizer MPS vs its generators, N = 4" begin
    rng = Xoshiro(41)
    for _ in 1:3
      m, S = stabilizer_mps(rng, 4, 25; maxdist = 3)
      v = mps_dense(m)
      for str in DR.all_paulis(4)
        @test DR.expect(v, str) ≈ SR.stab_expect(S, seam_strop(str)) atol = 1e-12
      end
    end
  end

  @testset "N = 46 Clifford states vs stabilizer reference" begin
    N = 46
    rng = Xoshiro(4646)
    for trial in 1:(SEAM_SLOW ? 4 : 1), d in (false, true)
      m, Smps = stabilizer_mps(rng, N, 120)
      @test maxlinkdim(m) >= 4                      # a real MPS, not a product state
      C = random_clifford(rng, N)
      Sψ = C * Smps                                  # generators of |ψ⟩ = C|mps⟩ (§4)
      # Q-space strings: high-weight random ones (⟨·⟩ mostly 0), the identity, prefix families,
      # and elements of the mps stabilizer group (⟨Q⟩ = ±1, so the frame sign matters).
      Qs = Set(designed_Qs(rng, N, 40))
      for _ in 1:20
        Q = prod(Smps[i] for i in rand(rng, 1:N, rand(rng, 1:6)))
        push!(Qs, seam_opstring(Q))
      end
      O = sum_from_Q(rng, C, collect(Qs); T = Float64)
      expected = sum(c * SR.stab_expect(Sψ, seam_strop(termstring(t, N))) for (t, c) in O.terms)
      @test count(!iszero, (SR.stab_expect(Sψ, seam_strop(termstring(t, N))) for t in keys(O.terms))) >= 10
      ψ = camps_from(m, C; daggered = d)
      @test abs(mt.expectation(ψ, O) - expected) <= EXACT_ATOL
      if SEAM_ITEM4
        evp = mt.expectation(ψ, O; method = :perterm)
        evm = mt.expectation(ψ, O; method = :mpo)
        @test abs(evp - expected) <= EXACT_ATOL
        @test abs(evp - evm) <= 1e-12
      else
        item4_pending()
      end
    end
  end

  # -- [item 4] keywords, empty sum, ψ unchanged, repeat calls ----------------------------

  @testset "[item 4] method, normalized, empty sum vs dense, N = $N" for N in (2, 3, 5, 8)
    if !SEAM_ITEM4
      item4_pending()
    else
      rng = Xoshiro(hash((:item4kw, N)))
      for trial in 1:ntrials, d in (false, true)
        m0 = random_raw_mps(rng, N, rand(rng, 1:8))
        v = mps_dense(m0)
        norm2 = real(sum(abs2, v))
        C = random_clifford(rng, N)
        for g in (:raw, (:centre, rand(rng, 1:N)))
          ψ = camps_from(in_gauge(m0, g), C; daggered = d)
          for T in (Float64, ComplexF64)
            O = sum_from_Q(rng, C, designed_Qs(rng, N, min(40, 4^N)); T)
            spec = spec_seam(v, inv(C), O)
            ev_default = mt.expectation(ψ, O)
            ev = Dict(meth => mt.expectation(ψ, O; method = meth) for meth in (:perterm, :mpo))
            @test isequal(ev_default, ev[:perterm])                # :perterm is the default
            @test abs(ev[:perterm] - spec) <= EXACT_ATOL
            @test abs(ev[:mpo] - spec) <= EXACT_ATOL
            @test abs(ev[:perterm] - ev[:mpo]) <= 1e-12
            for meth in (:perterm, :mpo)
              @test isequal(mt.expectation(ψ, O; method = meth, normalized = false), ev[meth])
              evn = mt.expectation(ψ, O; method = meth, normalized = true)
              @test evn isa ComplexF64
              @test abs(evn - spec / norm2) <= EXACT_ATOL
            end
            @test abs(norm2 - 1) > 1e-3                            # normalisation is visible
          end
          # Empty sum: 0, for both methods and normalisations.
          E = pp.PauliSum(Float64, N)
          for meth in (:perterm, :mpo), nz in (false, true)
            @test mt.expectation(ψ, E; method = meth, normalized = nz) == 0
          end
          @test mt.expectation(ψ, E) == 0
        end
      end
    end
  end

  @testset "[item 4] unknown method is an ArgumentError" begin
    if !SEAM_ITEM4
      item4_pending()
    else
      ψ = camps_from(random_raw_mps(Xoshiro(7), 3, 2), random_clifford(Xoshiro(8), 3); daggered = false)
      O = seam_ppsum(3, [("XYZ", 1.0)])
      for meth in (:dense, :MPO, :per_term, :exact)
        @test_throws ArgumentError mt.expectation(ψ, O; method = meth)
      end
    end
  end

  @testset "[item 4] ψ unchanged and repeat calls bitwise equal, N = $N" for N in (2, 4, 7)
    if !SEAM_ITEM4
      item4_pending()
    else
      rng = Xoshiro(hash((:item4unchanged, N)))
      for d in (false, true), g in (:raw, (:centre, 1), (:centre, N), (:centre, rand(rng, 1:N)))
        m0 = random_raw_mps(rng, N, rand(rng, 2:8))
        C = random_clifford(rng, N)
        ψ = camps_from(in_gauge(m0, g), C; daggered = d)
        snap = snapshot(ψ)
        O = sum_from_Q(rng, C, designed_Qs(rng, N, min(30, 4^N)); T = ComplexF64)
        first_ev = mt.expectation(ψ, O)
        @test unchanged(ψ, snap)
        for meth in (:perterm, :mpo), nz in (false, true)
          ev1 = mt.expectation(ψ, O; method = meth, normalized = nz)
          @test unchanged(ψ, snap)
          ev2 = mt.expectation(ψ, O; method = meth, normalized = nz)
          @test isequal(ev1, ev2)
        end
        @test isequal(mt.expectation(ψ, O), first_ev)
        redirect_stdout(devnull) do
          mt.expectation(ψ, O; verbose = true)
          mt.expectation(ψ, O; method = :mpo, verbose = true)
        end
        @test unchanged(ψ, snap)
        # Negative check of the helper: changing one tensor entry, or flipping the storage,
        # is noticed.
        ψbad = camps_from(in_gauge(m0, g), C; daggered = d)
        A = ψbad.mps[1]
        ψbad.mps[1] = A + 1e-15 * ITensor(randn(rng, ComplexF64, dim.(inds(A))...), inds(A)...)
        @test !unchanged(ψbad, snap)
        ψflip = camps_from(in_gauge(m0, g), C; daggered = !d)
        @test !unchanged(ψflip, snap)
      end
    end
  end

  @testset "[item 4] verbose = true returns the same value" begin
    if !SEAM_ITEM4
      item4_pending()
    else
      rng = Xoshiro(55)
      C = random_clifford(rng, 4)
      ψ = camps_from(random_raw_mps(rng, 4, 4), C; daggered = true)
      O = sum_from_Q(rng, C, designed_Qs(rng, 4, 30); T = ComplexF64)
      for meth in (:perterm, :mpo), nz in (false, true)
        ev = mt.expectation(ψ, O; method = meth, normalized = nz)
        evv = redirect_stdout(devnull) do
          mt.expectation(ψ, O; method = meth, normalized = nz, verbose = true)
        end
        @test isequal(evv, ev)
      end
    end
  end
end
