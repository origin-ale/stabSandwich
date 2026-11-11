using DisentangleCAMPS: apply, paulinature, disentangler, reducetoX, build_D, findfirstfreeXY, swapqubits
using QuantumClifford
using CliffordMPS
using Random
using Test

# Expected values follow test/CONVENTIONS.md §4 and are derived by hand, not from the code.
# Products use `A * B` = "B acts first" and `C * P` = C P C† (pinned in TestGoldenSet.jl).

# Controlled gates, first qubit = control. Hand-derived tableaus (X_c, X_t, Z_c, Z_t):
#   CY: X_c → X_c Y_t,  X_t → Z_c X_t,  Z_c → Z_c,  Z_t → Z_c Z_t
#   CZ: X_c → X_c Z_t,  X_t → Z_c X_t,  Z_c → Z_c,  Z_t → Z_t
const tCX = tCNOT
const tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
const tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

"The gate `gate` acting on qubits `inds` of an N-qubit register."
function embed(gate, N, inds)
  D = one(CliffordOperator, N)
  apply!(D, gate, inds)
  return D
end

@testset "assumption pins: controlled gates and embedding" begin
  @test tCY == C"XY ZX ZI ZZ"
  @test tCZ == C"XZ ZX ZI IZ"
  # embed(gate, N, [c, t]): the first index is the control.
  @test embed(tCNOT, 2, [1, 2]) == tCNOT
  @test embed(tCNOT, 2, [2, 1]) * P"IX" == P"+XX"
  @test embed(tCNOT, 2, [2, 1]) * P"XI" == P"+XI"
  @test embed(tPhase, 3, [2]) * P"XXX" == P"+XYX"
  @test embed(tSWAP, 3, [1, 3]) * P"XYZ" == P"+ZYX"
end

@testset "apply!(clifford, op)" begin
  C0 = one(CliffordOperator, 3) # Identity

  C = deepcopy(C0)
  @assert C == C"XII IXI IIX ZII IZI IIZ" "Initialization or Clifford macro not working as expected!"
  apply!(C, tHadamard, [1])
  @test C == C"ZII IXI IIX XII IZI IIZ"
  apply!(C, tCNOT, [1,3])
  @test C == C"ZII IXI IIX XIX IZI ZIZ"
  apply!(C, tCNOT, [2,1])
  @test C == C"ZZI XXI IIX XIX IZI ZZZ"
  apply!(C, tPhase, [2])
  @test C == C"ZZI XYI IIX XIX IZI ZZZ"
  apply!(C, tHadamard, [2])
  @test C == C"ZXI -XYI IIX XIX IXI ZXZ"

  C = deepcopy(C0)
  @test apply!(C, tHadamard, [2]) === nothing
  @test C != C0

  C = deepcopy(C0)
  apply!(C, tHadamard, [2])
  apply!(C, tHadamard, [2])
  @test C == C0

  C = deepcopy(C0)
  apply!(C, tCNOT, [1, 3])
  apply!(C, tCNOT, [1, 3])
  @test C == C0

  C = deepcopy(C0)
  for _ in 1:4
    apply!(C, tPhase, [2])
  end
  @test C == C0

  op = one(CliffordOperator, 2)
  apply!(op, tHadamard, [1])
  apply!(op, tCNOT, [1, 2])

  C1 = deepcopy(C0)
  apply!(C1, op, [3, 1])

  C2 = deepcopy(C0)
  apply!(C2, tHadamard, [3])
  apply!(C2, tCNOT, [3, 1])

  @test C1 == C2

  @test_throws BoundsError apply!(deepcopy(C0), tHadamard, [0])
  @test_throws BoundsError apply!(deepcopy(C0), tHadamard, [1, 2])

  # Current behavior in QuantumClifford does not reject these malformed indices.
  @test_broken try
    apply!(deepcopy(C0), tHadamard, [4])
    false
  catch e
    e isa BoundsError
  end

  @test_broken try
    apply!(deepcopy(C0), tCNOT, [1])
    false
  catch e
    e isa BoundsError
  end
end

@testset "apply(pauli, cliff)" begin
  # apply(P, C) = C† P C (CONVENTIONS.md §4).
  Cid = one(CliffordOperator, 3)
  @test apply(P"YZI", Cid) == P"+YZI"

  # C = S:  S† X S = -Y,  S† Y S = +X,  S† Z S = Z.
  C1 = embed(tPhase, 1, [1])
  @test apply(P"X", C1) == P"-Y"
  @test apply(P"Y", C1) == P"+X"
  @test apply(P"Z", C1) == P"+Z"

  # C = H on qubit 1 (self-inverse).
  C2 = embed(tHadamard, 2, [1])
  @test apply(P"XZ", C2) == P"+ZZ"
  @test apply(P"ZZ", C2) == P"+XZ"

  # C = CNOT 1→2 (self-inverse).
  C3 = embed(tCNOT, 2, [1, 2])
  @test apply(P"XI", C3) == P"+XX"
  @test apply(P"IZ", C3) == P"+ZZ"
  @test apply(P"YI", C3) == P"+YX"
  @test apply(P"IY", C3) == P"+ZY"

  # C = S₁ · CNOT(3→1) · H₃ (H₃ acts first), so C† P C conjugates P by S₁†, then CNOT, then H₃.
  #   XIX → -YIX → -ZIY → +ZIY
  #   ZZI →  ZZI →  ZZZ → +ZZX
  #   IYZ →  IYZ →  IYZ → +IYX
  C4 = embed(tPhase, 3, [1]) * embed(tCNOT, 3, [3, 1]) * embed(tHadamard, 3, [3])
  @test apply(P"XIX", C4) == P"+ZIY"
  @test apply(P"ZZI", C4) == P"+ZZX"
  @test apply(P"IYZ", C4) == P"+IYX"

  # Not self-inverse, so the direction matters: C† P C ≠ C P C†.
  @test apply(P"XIX", C4) != C4 * P"XIX"
  @test apply(P"XIX", C4) == inv(C4) * P"XIX"

  @test_throws DimensionMismatch apply(P"XYZ", one(CliffordOperator, 2))
end

@testset "paulinature on C=I" begin
  dnt = :disentanglable
  log = :logical
  trv = :trivial

  C = one(CliffordOperator, 3)
  @test paulinature(0, C, P"IXI") == dnt
  @test paulinature(1, C, P"IXI") == dnt
  @test paulinature(2, C, P"IXI") == log

  @test paulinature(0, C, P"IYI") == dnt
  @test paulinature(1, C, P"IYI") == dnt
  @test paulinature(2, C, P"IYI") == log

  @test paulinature(0, C, P"IZI") == trv
  @test paulinature(1, C, P"IZI") == trv
  @test paulinature(2, C, P"IZI") == log

  @test paulinature(0, C, P"XXI") == dnt
  @test paulinature(1, C, P"XXI") == dnt
  @test paulinature(2, C, P"XXI") == log

  @test paulinature(0, C, P"XZI") == dnt
  @test paulinature(1, C, P"XZI") == log
  @test paulinature(2, C, P"XZI") == log

  @test paulinature(0, C, P"XIZ") == dnt
  @test paulinature(1, C, P"XIZ") == log
  @test paulinature(2, C, P"XIZ") == log
end

@testset "paulinature on C≠I" begin
  dnt = :disentanglable
  log = :logical
  trv = :trivial

  # C1 = H on qubit 1: XZI -> ZZI and ZXI -> XXI.
  C1 = one(CliffordOperator, 3)
  apply!(C1, tHadamard, [1])

  @test paulinature(0, C1, P"XZI") == trv
  @test paulinature(1, C1, P"XZI") == log
  @test paulinature(2, C1, P"XZI") == log

  @test paulinature(0, C1, P"ZXI") == dnt
  @test paulinature(1, C1, P"ZXI") == dnt
  @test paulinature(2, C1, P"ZXI") == log

  # C2 = S on qubits 1 and 3: XIX -> YIY and ZZZ -> ZZZ.
  C2 = one(CliffordOperator, 3)
  apply!(C2, tPhase, [1])
  apply!(C2, tPhase, [3])

  @test paulinature(0, C2, P"XIX") == dnt
  @test paulinature(1, C2, P"XIX") == dnt
  @test paulinature(2, C2, P"XIX") == dnt

  @test paulinature(0, C2, P"ZZZ") == trv
  @test paulinature(1, C2, P"ZZZ") == log
  @test paulinature(2, C2, P"ZZZ") == log

  # C3 = CNOT 1->2: XIZ -> XXZ and IZI -> ZZI.
  C3 = one(CliffordOperator, 3)
  apply!(C3, tCNOT, [1, 2])

  @test paulinature(0, C3, P"XIZ") == dnt
  @test paulinature(1, C3, P"XIZ") == dnt
  @test paulinature(2, C3, P"XIZ") == log

  @test paulinature(0, C3, P"IZI") == trv
  @test paulinature(1, C3, P"IZI") == log
  @test paulinature(2, C3, P"IZI") == log

  # Edge case k = n: upper region is empty.
  Cid = one(CliffordOperator, 3)
  @test paulinature(3, Cid, P"III") == trv
  @test paulinature(3, Cid, P"ZZZ") == log
end

@testset "swapqubits" begin
  C = one(CliffordOperator, 3)

  @test swapqubits(P"XYZ", 1, 3)[1] == P"ZYX"

  @test swapqubits(P"IIX", 2, 3)[1] == P"IXI"
end

@testset "findfirstfreeXY" begin
  # X on a free qubit (k = number of magic qubits at the start)
  @test findfirstfreeXY(P"IXI", 1) == 2
  @test findfirstfreeXY(P"IIX", 2) == 3

  # Y on a free qubit (Y has xbit=true in the Pauli representation)
  @test findfirstfreeXY(P"IYI", 1) == 2
  @test findfirstfreeXY(P"IIY", 2) == 3

  # X on multiple free qubits: return the first one
  @test findfirstfreeXY(P"IXX", 1) == 2

  # k = 0: scan from qubit 1
  @test findfirstfreeXY(P"XIZI", 0) == 1
  @test findfirstfreeXY(P"ZIXI", 2) == 3

  # Single free qubit
  @test findfirstfreeXY(P"IXI", 1) == 2
end

@testset "reducetoX" begin
  # Spec: the gate is the phase gate S on qubit i, which turns Y_i into X_i.
  # The sign of the returned string is bookkeeping and is not pinned
  # (CONVENTIONS.md §4); only its X/Z content is.
  samebits(P, Q) = xbit(P) == xbit(Q) && zbit(P) == zbit(Q)

  reduced, phase = reducetoX(P"XYZ", 2)
  @test samebits(reduced, P"XXZ")
  @test phase == embed(tPhase, 3, [2])

  reduced, phase = reducetoX(P"YYY", 1)
  @test samebits(reduced, P"XYY")
  @test phase == embed(tPhase, 3, [1])
end

@testset "build_D" begin
  test_D = one(CliffordOperator, 5)
  apply!(test_D, tCY, [3,1])
  apply!(test_D, tCX, [3,2])
  apply!(test_D, tCZ, [3,4])
  apply!(test_D, tCX, [3,5])
  @test build_D(P"YXXZX", 3) == test_D
end

# CONVENTIONS.md §4, correctness properties of (D, s) for Q = C† P C and k magic qubits.
"Property 1: D† Q D = s · X_{k+1}."
function maps_to_X(D, s, Q, k)
  N = length(Q)
  target = s == 1 ? single_x(N, k + 1) : -single_x(N, k + 1)
  return abs(s) == 1 && inv(D) * Q == target
end

"Is g a +1-signed product of Z's on qubits k+1…N only?"
function in_free_Z_group(g, k)
  return g.phase[] == 0x0 && !any(xbit(g)) && !any(zbit(g)[1:k])
end

"""Property 2: D fixes every state |φ⟩₁…ₖ ⊗ |0…0⟩ (up to a global phase).
Equivalent tableau conditions, with G the group generated by Z_{k+1}, …, Z_N:
D Z_j D† ∈ G for j > k, and D A D† ∈ A·G for A = X_j, Z_j with j ≤ k."""
function fixes_magic_register(D, k)
  N = length(D[1])
  ok = true
  for j in 1:k, A in (single_x(N, j), single_z(N, j))
    ok &= in_free_Z_group((D * A) * A, k)
  end
  for j in k+1:N
    ok &= in_free_Z_group(D * single_z(N, j), k)
  end
  return ok
end

@testset "property helpers reject wrong disentanglers" begin
  # Q = YXI, k = 1; the valid D is CY(2→1) with s = +1.
  Q = P"YXI"
  D = embed(tCY, 3, [2, 1])
  @test maps_to_X(D, +1, Q, 1) && fixes_magic_register(D, 1)
  @test !maps_to_X(D, -1, Q, 1)                              # wrong sign
  @test !maps_to_X(embed(tCX, 3, [2, 1]), +1, Q, 1)          # wrong controlled Pauli
  @test !maps_to_X(one(CliffordOperator, 3), +1, Q, 1)       # no disentangling at all
  # Property 2 on its own:
  @test !fixes_magic_register(embed(tHadamard, 3, [3]), 1)   # rotates a free |0⟩
  @test !fixes_magic_register(embed(tHadamard, 3, [1]), 1)   # acts on the magic qubit
  @test !fixes_magic_register(embed(tCX, 3, [1, 2]), 1)      # control on the magic qubit
  @test !fixes_magic_register(embed(tSWAP, 3, [1, 2]), 1)    # swaps magic and free
  @test fixes_magic_register(embed(tSWAP, 3, [2, 3]), 1)     # swaps two free |0⟩: allowed
  @test fixes_magic_register(embed(tPhase, 3, [3]), 1)       # S on a free |0⟩: allowed
end

@testset "disentangler: pinned construction" begin
  # Construction (CONVENTIONS.md §4): i = smallest index > k with X or Y in Q;
  # D₁ = S_i if Q_i = Y;  D₂ = SWAP(i, k+1) if i ≠ k+1;
  # D₃ = controlled-Q′_j from k+1 onto every other non-identity qubit j of
  # Q′ = (D₁D₂)† Q (D₁D₂);  D = D₁ D₂ D₃;  s = sign of Q.
  function check(k, C, P, Q_hand, D_hand, s_hand)
    Q = inv(C) * P
    @test Q == Q_hand                       # the hand-derived Q itself
    @test paulinature(k, C, P) == :disentanglable
    D, s = disentangler(k, C, P)
    @test D == D_hand
    @test s == s_hand
    @test maps_to_X(D_hand, s_hand, Q, k)   # the hand-derived D is a valid disentangler
    @test fixes_magic_register(D_hand, k)
  end

  SW12 = C"IXI XII IIX IZI ZII IIZ"         # SWAP(1,2), self-inverse
  @assert SW12 == embed(tSWAP, 3, [1, 2])

  # Q = YXI, k = 1: i = 2 = k+1, X there. D₃ = CY(2→1).
  check(1, SW12, P"XYI", P"+YXI", embed(tCY, 3, [2, 1]), +1)

  # Q = YYI, k = 1: i = 2, Y there → D₁ = S₂, Q′ = YXI. D₃ = CY(2→1).
  check(1, SW12, P"YYI", P"+YYI", embed(tPhase, 3, [2]) * embed(tCY, 3, [2, 1]), +1)

  # Q = YIX, k = 1: i = 3 ≠ k+1 → D₂ = SWAP(2,3), Q′ = YXI. D₃ = CY(2→1).
  check(1, SW12, P"IYX", P"+YIX", embed(tSWAP, 3, [2, 3]) * embed(tCY, 3, [2, 1]), +1)

  # C = 1, k = 0, Q = ZYX: i = 2, Y → D₁ = S₂ (ZXX); D₂ = SWAP(1,2) (XZX);
  # D₃ = CZ(1→2) · CX(1→3).
  Id3 = one(CliffordOperator, 3)
  check(0, Id3, P"ZYX", P"+ZYX",
        embed(tPhase, 3, [2]) * embed(tSWAP, 3, [1, 2]) * embed(tCZ, 3, [1, 2]) * embed(tCX, 3, [1, 3]), +1)

  # Negative sign: C = 1, k = 0, Q = -XZ: i = 1, X. D₃ = CZ(1→2), s = -1.
  check(0, one(CliffordOperator, 2), P"-XZ", P"-XZ", embed(tCZ, 2, [1, 2]), -1)

  # C = H₃ (self-inverse), k = 2, N = 4. P = ZYZZ → Q = ZYXZ: i = 3 = k+1, X.
  # D₃ = CZ(3→1) · CY(3→2) · CZ(3→4).
  H3 = embed(tHadamard, 4, [3])
  check(2, H3, P"ZYZZ", P"+ZYXZ",
        embed(tCZ, 4, [3, 1]) * embed(tCY, 4, [3, 2]) * embed(tCZ, 4, [3, 4]), +1)

  # Sign from the frame: C = H₃, k = 2, P = IIYI → Q = -IIYI (H Y H = -Y).
  # i = 3, Y → D₁ = S₃; nothing else to do. s = -1.
  check(2, H3, P"IIYI", P"-IIYI", embed(tPhase, 4, [3]), -1)

  # Y at i, a swap and magic-qubit targets together: C = 1, k = 1, N = 4, Q = -XIZY.
  # i = 4, Y → D₁ = S₄ (XIZX); D₂ = SWAP(4,2) (XXZI); D₃ = CX(2→1) · CZ(2→3). s = -1.
  check(1, one(CliffordOperator, 4), P"-XIZY", P"-XIZY",
        embed(tPhase, 4, [4]) * embed(tSWAP, 4, [4, 2]) * embed(tCX, 4, [2, 1]) * embed(tCZ, 4, [2, 3]), -1)
end

@testset "disentangler: correctness properties (seeded random)" begin
  rng = Xoshiro(20261111)
  tested = 0
  for N in 2:6, _ in 1:60
    k = rand(rng, 0:N-1)
    C = random_clifford(rng, N)
    P = random_pauli(rng, N; nophase = false, realphase = true)
    Q = inv(C) * P
    any(xbit(Q)[k+1:N]) || continue          # only the disentanglable case
    tested += 1
    D, s = disentangler(k, C, P)
    @test s == (Q.phase[] == 0x0 ? 1 : -1)
    @test maps_to_X(D, s, Q, k)
    @test fixes_magic_register(D, k)
    @test paulinature(k, C, P) == :disentanglable
  end
  @test tested >= 200                         # guard against a vacuous pass
end
