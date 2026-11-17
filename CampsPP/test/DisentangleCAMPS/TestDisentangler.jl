using CampsPP.DisentangleCAMPS: _apply_on_sites!, apply, paulinature, disentangler, reducetoX, build_D, findfirstfreeXY, swapqubits
using QuantumClifford
using Random
using Test

# Expected values follow test/CONVENTIONS.md §4 and are derived by hand, not from the code.
# Products use `A * B` = "B acts first" and `cliffconj(C, P)` = C P C† (helpers.jl; both
# pinned in TestGoldenSet.jl).

# Controlled gates, first qubit = control. Hand-derived tableaus (X_c, X_t, Z_c, Z_t):
#   CY: X_c → X_c Y_t,  X_t → Z_c X_t,  Z_c → Z_c,  Z_t → Z_c Z_t
#   CZ: X_c → X_c Z_t,  X_t → Z_c X_t,  Z_c → Z_c,  Z_t → Z_t
const tCX = tCNOT
const tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
const tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

@testset "assumption pins: controlled gates and embedding" begin
  @test tCY == C"XY ZX ZI ZZ"
  @test tCZ == C"XZ ZX ZI IZ"
  # embed(gate, N, [c, t]): the first index is the control.
  @test embed(tCNOT, 2, [1, 2]) == tCNOT
  @test cliffconj(embed(tCNOT, 2, [2, 1]), P"IX") == P"+XX"
  @test cliffconj(embed(tCNOT, 2, [2, 1]), P"XI") == P"+XI"
  @test cliffconj(embed(tPhase, 3, [2]), P"XXX") == P"+XYX"
  @test cliffconj(embed(tSWAP, 3, [1, 3]), P"XYZ") == P"+ZYX"
  @test cliffconj(embed(tPhase, 3, [2]), P"XYX") == P"-XXX"   # S Y S† = −X: signs carried
  @test cliffconj(embed(tCY, 3, [3, 1]), P"XII") == P"+XIZ"     # X_t → Z_c X_t (control 3)
  @test cliffconj(embed(tCY, 3, [3, 1]), P"IIX") == P"+YIX"     # X_c → Y_t X_c
  @test_throws ArgumentError embed(tCNOT, 3, [2, 2])
  @test_throws ArgumentError embed(tCNOT, 3, [1, 4])
end

# `_apply_on_sites!(C, op, sites)` (C ← op·C), the private replacement for the vendored
# `apply!(::CliffordOperator, op, sites)` that swapqubits, reducetoX and build_D use (plan,
# design correction 3). Expected values unchanged from the old tests of the vendored method,
# plus a seeded check against the QC-native `embed` (helpers.jl).
@testset "_apply_on_sites!(clifford, op, sites)" begin
  C0 = one(CliffordOperator, 3) # Identity

  C = deepcopy(C0)
  @assert C == C"XII IXI IIX ZII IZI IIZ" "Initialization or Clifford macro not working as expected!"
  _apply_on_sites!(C, tHadamard, [1])
  @test C == C"ZII IXI IIX XII IZI IIZ"
  _apply_on_sites!(C, tCNOT, [1,3])
  @test C == C"ZII IXI IIX XIX IZI ZIZ"
  _apply_on_sites!(C, tCNOT, [2,1])
  @test C == C"ZZI XXI IIX XIX IZI ZZZ"
  _apply_on_sites!(C, tPhase, [2])
  @test C == C"ZZI XYI IIX XIX IZI ZZZ"
  _apply_on_sites!(C, tHadamard, [2])
  @test C == C"ZXI -XYI IIX XIX IXI ZXZ"

  C = deepcopy(C0)
  @test _apply_on_sites!(C, tHadamard, [2]) === C          # in place, returns C
  @test C != C0

  C = deepcopy(C0)
  _apply_on_sites!(C, tHadamard, [2])
  _apply_on_sites!(C, tHadamard, [2])
  @test C == C0

  C = deepcopy(C0)
  _apply_on_sites!(C, tCNOT, [1, 3])
  _apply_on_sites!(C, tCNOT, [1, 3])
  @test C == C0

  C = deepcopy(C0)
  for _ in 1:4
    _apply_on_sites!(C, tPhase, [2])
  end
  @test C == C0

  op = one(CliffordOperator, 2)
  _apply_on_sites!(op, tHadamard, [1])
  _apply_on_sites!(op, tCNOT, [1, 2])

  C1 = deepcopy(C0)
  _apply_on_sites!(C1, op, [3, 1])

  C2 = deepcopy(C0)
  _apply_on_sites!(C2, tHadamard, [3])
  _apply_on_sites!(C2, tCNOT, [3, 1])

  @test C1 == C2

  # C ← op·C with op placed on `sites`, for random C and op (1- and 2-qubit, any order).
  rng = Xoshiro(20261117)
  for N in 2:5, _ in 1:20
    n = rand(rng, 1:2)
    C, op, sites = random_clifford(rng, N), random_clifford(rng, n), randperm(rng, N)[1:n]
    @test _apply_on_sites!(copy(C), op, sites) == embed(op, N, sites) * C
  end

  @test_throws BoundsError _apply_on_sites!(deepcopy(C0), tHadamard, [0])
  @test_throws BoundsError _apply_on_sites!(deepcopy(C0), tHadamard, [1, 2])

  # Not checked: an index beyond N ([4]) or too few indices ([1] for a 2-qubit gate).
  # QuantumClifford's `apply!(::Stabilizer, op, sites)` does not reject them, and the callers
  # only pass indices built from a valid Q (the 2 `@test_broken` that pinned this in the
  # vendored method were dropped in Phase 6).
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
  @test apply(P"XIX", C4) != cliffconj(C4, P"XIX")
  @test apply(P"XIX", C4) == cliffconj(inv(C4), P"XIX")

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
  C1 = embed(tHadamard, 3, [1])

  @test paulinature(0, C1, P"XZI") == trv
  @test paulinature(1, C1, P"XZI") == log
  @test paulinature(2, C1, P"XZI") == log

  @test paulinature(0, C1, P"ZXI") == dnt
  @test paulinature(1, C1, P"ZXI") == dnt
  @test paulinature(2, C1, P"ZXI") == log

  # C2 = S on qubits 1 and 3: XIX -> YIY and ZZZ -> ZZZ.
  C2 = embed(tPhase, 3, [3]) * embed(tPhase, 3, [1])

  @test paulinature(0, C2, P"XIX") == dnt
  @test paulinature(1, C2, P"XIX") == dnt
  @test paulinature(2, C2, P"XIX") == dnt

  @test paulinature(0, C2, P"ZZZ") == trv
  @test paulinature(1, C2, P"ZZZ") == log
  @test paulinature(2, C2, P"ZZZ") == log

  # C3 = CNOT 1->2: XIZ -> XXZ and IZI -> ZZI.
  C3 = embed(tCNOT, 3, [1, 2])

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

@testset "paulinature: direction-sensitive frame" begin
  # C4 = S₁ · CNOT(3→1) · H₃ is not self-inverse, so C†PC and CPC† differ. By hand
  # (see "apply(pauli, cliff)"): Q(XIX) = +ZIY, Q(ZZI) = +ZZX, Q(IYZ) = +IYX.
  # And Q(YIX) = IIZ:  YIX → (S₁†) XIX → (CNOT 3→1) IIX → (H₃) IIZ.
  # With CPC† instead: XIX ↦ YIZ, which would be :logical already at k = 1.
  dnt, log, trv = :disentanglable, :logical, :trivial
  C4 = embed(tPhase, 3, [1]) * embed(tCNOT, 3, [3, 1]) * embed(tHadamard, 3, [3])
  @assert cliffconj(inv(C4), P"YIX") == P"+IIZ"

  for P in (P"XIX", P"ZZI", P"IYZ")              # X or Y on qubit 3 only
    @test [paulinature(k, C4, P) for k in 0:3] == [dnt, dnt, dnt, log]
  end
  @test [paulinature(k, C4, P"YIX") for k in 0:3] == [trv, trv, trv, log]   # Q = IIZ
  @test paulinature(2, C4, P"-XIX") == dnt       # the sign of P does not matter
  @test paulinature(3, C4, P"-YIX") == log
end

@testset "paulinature: seeded random against the spec rule" begin
  rng = Xoshiro(20261112)
  seen = Dict(:disentanglable => 0, :logical => 0, :trivial => 0)
  mismatches = 0
  for N in 2:6, _ in 1:80
    k = rand(rng, 0:N)
    C = random_clifford(rng, N)
    # Half of the draws are built to be non-disentanglable: Q is I/Z beyond qubit k.
    Q = random_pauli(rng, N; nophase = false, realphase = true)
    if rand(rng, Bool)
      for j in k+1:N
        Q[j] = (false, Q[j][2])
      end
    end
    P = cliffconj(C, Q)                                   # so that C†PC = Q
    nature = spec_nature(k, Q)
    seen[nature] += 1
    mismatches += paulinature(k, C, P) != nature
  end
  @test mismatches == 0
  @test all(>=(40), values(seen))                # every branch is well represented
end

@testset "swapqubits" begin
  C = one(CliffordOperator, 3)

  @test swapqubits(P"XYZ", 1, 3)[1] == P"ZYX"

  @test swapqubits(P"IIX", 2, 3)[1] == P"IXI"

  # The returned Clifford is the SWAP itself.
  @test swapqubits(P"XYZ", 1, 3)[2] == embed(tSWAP, 3, [1, 3])
  @test swapqubits(P"IIX", 2, 3)[2] == embed(tSWAP, 3, [2, 3])
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
  # CY(3→1), then CX(3→2), CZ(3→4), CX(3→5) (they commute: shared control).
  test_D = embed(tCX, 5, [3, 5]) * embed(tCZ, 5, [3, 4]) * embed(tCX, 5, [3, 2]) * embed(tCY, 5, [3, 1])
  @test build_D(P"YXXZX", 3) == test_D
end

# CONVENTIONS.md §4, correctness properties of (D, s) for Q = C† P C and k magic qubits.
"Property 1: D† Q D = s · X_{k+1}."
function maps_to_X(D, s, Q, k)
  N = length(Q)
  target = s == 1 ? single_x(N, k + 1) : -single_x(N, k + 1)
  return abs(s) == 1 && cliffconj(inv(D), Q) == target
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
    ok &= in_free_Z_group(cliffconj(D, A) * A, k)
  end
  for j in k+1:N
    ok &= in_free_Z_group(cliffconj(D, single_z(N, j)), k)
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
    Q = cliffconj(inv(C), P)
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
    Q = cliffconj(inv(C), P)
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
