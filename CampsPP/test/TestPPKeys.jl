# Item 3 (post-port plan): PauliPropagation term integers (keys) — assumption pins for the
# PP ≥ 0.9 upgrade, written before it. Expected values come from test/CONVENTIONS.md §3 and
# hand values only.
#
# - "key type rule (PP ≥ 0.9)" pins `getinttype` to the §3 rule (it failed on 0.7.3 by
#   design: 0.7.3 used the 8⌈2N/8⌉-bit types).
# - The bit layout and the integer operations that the converter (`pauliop_from_ppterm`, via
#   `pp.getpauli`) and the tests (`termstring` in helpers.jl, TestInitialize.jl) use on keys
#   are the same on both versions, so they run on any version. They cover every type boundary
#   of the ≥ 0.9 rule (N = 4|5, 8|9, 16|17, 32|33, 64|65, 96|97, 128|129) and the 64-bit limb
#   boundaries above 64 qubits (qubits 32|33, 64|65, 96|97, 128|129).

# -- Spec (CONVENTIONS.md §3) -------------------------------------------------------------

"Code of a Pauli character: I, X, Y, Z = 0, 1, 2, 3."
paulicode(c::Char) = findfirst(==(c), "IXYZ") - 1

"The §3 term of an N-character string (character j = qubit j), as a `BigInt`: Σⱼ codeⱼ · 4^(j−1)."
spec_term(str::AbstractString) = sum(big(paulicode(c)) << (2 * (j - 1)) for (j, c) in enumerate(str); init = big(0))

"""
The §3 key type at N qubits on PP ≥ 0.9: the smallest of UInt8 … UInt128 that holds 2N bits,
above that `NTupleInteger{⌈2N/64⌉}`.
"""
function spec_keytype(N)
  for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
    2N <= 8 * sizeof(T) && return T
  end
  return pp.NTupleInteger{cld(2N, 64)}
end

"A pp string from an N-character string, every qubit named (I included), coefficient c."
ppstring(str, c = 1.0) = pp.PauliString(length(str), Symbol.(collect(str)), collect(1:length(str)), c)

# Type boundaries of the ≥ 0.9 rule, both sides; N = 10 (old UInt24) and 46 (production).
const BOUNDARY_N = [1, 4, 5, 8, 9, 10, 16, 17, 32, 33, 46, 64, 65, 70, 96, 97, 128, 129]

"Qubits next to a 64-bit word boundary (32|33, 64|65, …) that exist at N."
boundary_qubits(N) = filter(j -> 1 <= j <= N, vcat([[32w, 32w + 1] for w in 1:fld(N, 32)]...))

"Random N-character string with non-identity Paulis on every qubit next to a word boundary."
function random_boundary_string(rng, N)
  chars = [rand(rng) < 0.3 ? rand(rng, "XYZ") : 'I' for _ in 1:N]
  for j in boundary_qubits(N)
    chars[j] = rand(rng, "XYZ")
  end
  chars[N] = rand(rng, "XYZ")
  return String(chars)
end

@testset "Item 3: PauliPropagation keys (CONVENTIONS.md §3)" begin

  @testset "spec helpers (negative checks)" begin
    @test spec_term("I") == 0 && spec_term("X") == 1 && spec_term("Y") == 2 && spec_term("Z") == 3
    @test spec_term("IX") == 4 && spec_term("ZY") == 3 + 2 * 4      # qubit 1 lowest
    @test spec_term("XZ") != spec_term("ZX")                         # order matters
    @test spec_term("X") != spec_term("Z")                           # X and Z differ
    @test spec_term("I"^32 * "Z") == big(3) << 64                    # first qubit of word 2
    # The key-type helper against hand values of the rule (the 0.7.3 rule would give UInt24 at 10).
    @test [spec_keytype(N) for N in (4, 5, 8, 9, 10, 16, 17, 32, 33, 46, 64)] ==
          [UInt8, UInt16, UInt16, UInt32, UInt32, UInt32, UInt64, UInt64, UInt128, UInt128, UInt128]
    @test spec_keytype(10) != UInt16 && spec_keytype(46) != UInt64
    @test boundary_qubits(70) == [32, 33, 64, 65] && boundary_qubits(32) == [32]
  end

  @testset "key type rule (PP ≥ 0.9)" begin
    NT = pp.NTupleInteger
    # Hand values at every boundary of the rule (§3).
    hand = [1 => UInt8, 4 => UInt8, 5 => UInt16, 8 => UInt16, 9 => UInt32, 10 => UInt32,
            16 => UInt32, 17 => UInt64, 32 => UInt64, 33 => UInt128, 46 => UInt128,
            64 => UInt128, 65 => NT{3}, 70 => NT{3}, 96 => NT{3}, 97 => NT{4}, 128 => NT{4},
            129 => NT{5}]
    @testset "N = $N" for (N, T) in hand
      @test pp.getinttype(N) === T
      @test spec_keytype(N) === T
      @test keytype(pp.PauliSum(N).terms) === T
      @test ppstring("I"^(N - 1) * "Z").term isa T
    end
    @test all(pp.getinttype(N) === spec_keytype(N) for N in 1:200)
    @test all(pp.getinttype(N) <: Unsigned for N in 1:200)
    @test pp.getinttype(46) === UInt128                  # production size
    # `getinttype` never returns a 2-limb NTupleInteger (≤ 128 bits is UInt128).
    @test !any(pp.getinttype(N) === NT{2} for N in 1:200)
  end

  @testset "bit layout, N = $N" for N in BOUNDARY_N
    T = pp.getinttype(N)
    # Every single-qubit string, by hand from §3: code << 2(j−1).
    @test all(BigInt(ppstring("I"^(j - 1) * s * "I"^(N - j)).term) == big(paulicode(s)) << (2 * (j - 1))
              for j in 1:N, s in ('X', 'Y', 'Z'))
    @test iszero(ppstring("I"^N).term)
    # Random multi-qubit strings, with Paulis on both sides of each word boundary.
    rng = Xoshiro(hash((:item3layout, N)))
    @test all(1:50) do _
      str = random_boundary_string(rng, N)
      BigInt(ppstring(str).term) == spec_term(str)
    end
  end

  @testset "bit layout across 64-bit words: hand pins" begin
    # X on qubit 32 = bits 62–63 (code 1 → bit 62); Z on 33 = bits 64–65 (code 3).
    @test BigInt(pp.PauliString(70, [:X, :Z], [32, 33]).term) == (big(1) << 62) | (big(3) << 64)
    @test BigInt(pp.PauliString(46, [:X, :Z], [32, 33]).term) == (big(1) << 62) | (big(3) << 64)
    # Y on 64 = bits 126–127 (code 2 → bit 127); X on 65 = bit 128.
    @test BigInt(pp.PauliString(70, [:Y, :X], [64, 65]).term) == (big(2) << 126) | (big(1) << 128)
    @test BigInt(pp.PauliString(129, [:Z], [129]).term) == big(3) << 256
  end

  @testset "integer operations on keys, N = $N" for N in BOUNDARY_N
    T = pp.getinttype(N)
    rng = Xoshiro(hash((:item3ops, N)))
    str = random_boundary_string(rng, N)
    t = ppstring(str).term
    @test t isa Integer                                    # `pauliop_from_ppterm(term::Integer, …)`
    @test t isa Unsigned
    # `termstring` (helpers.jl): Int((t >> 2(j−1)) & 0x3) is the code of qubit j.
    @test [Int((t >> (2 * (j - 1))) & 0x3) for j in 1:N] == paulicode.(collect(str))
    @test termstring(t, N) == str
    # `pauliop_from_ppterm`: Int(pp.getpauli(t, j)) is the code of qubit j.
    @test [Int(pp.getpauli(t, j)) for j in 1:N] == paulicode.(collect(str))
    # TestInitialize.jl: T(c) << 2(j−1) is the key of c on qubit j, of type T, and finds
    # PP's own key in a Dict (equal and same hash).
    for j in unique([1; boundary_qubits(N); N]), (s, c) in ((:X, 1), (:Y, 2), (:Z, 3))
      k = T(c) << (2 * (j - 1))
      @test k isa T
      @test k == pp.PauliString(N, [s], [j]).term
      @test haskey(pp.PauliSum(pp.PauliString(N, [s], [j])).terms, k)
    end
    @test zero(T) == ppstring("I"^N).term
    @test haskey(pp.PauliSum(ppstring("I"^N, 0.5)).terms, zero(T))
    # Negative check: a shift by the wrong amount reads another qubit's code.
    j = N
    @test Int((t >> (2 * (j - 1))) & 0x3) != 0             # qubit N is not I (random_boundary_string)
    N > 1 && @test Int((t >> (2 * j)) & 0x3) == 0            # nothing above qubit N
  end
end
