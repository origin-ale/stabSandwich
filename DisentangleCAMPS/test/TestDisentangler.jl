using DisentangleCAMPS: apply, paulinature, disentangler, reducetoX, build_D, findfirstfreeXY, swapqubits
using QuantumClifford
using CliffordMPS
using Test

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
  C = C"-XX iZY -XZ IZ"
  @test apply(P"YY", C) == P"+ZZ"

  Cid = one(CliffordOperator, 3)
  @test apply(P"YZI", Cid) == P"+YZI"

  C1 = one(CliffordOperator, 1)
  apply!(C1, tPhase, [1])
  @test apply(P"X", C1) == P"+Y"
  @test apply(P"Y", C1) == P"-X"
  @test apply(P"Z", C1) == P"+Z"

  C2 = one(CliffordOperator, 2)
  apply!(C2, tHadamard, [1])
  @test apply(P"XZ", C2) == P"+ZZ"
  @test apply(P"ZZ", C2) == P"+XZ"

  C3 = one(CliffordOperator, 2)
  apply!(C3, tCNOT, [1, 2])
  @test apply(P"XI", C3) == P"+XX"
  @test apply(P"IZ", C3) == P"+ZZ"
  @test apply(P"YI", C3) == P"+YX"
  @test apply(P"IY", C3) == P"+ZY"

  C4 = one(CliffordOperator, 3)
  apply!(C4, tHadamard, [3])
  apply!(C4, tCNOT, [3, 1])
  apply!(C4, tPhase, [1])
  @test apply(P"XIX", C4) == P"+YIZ"
  @test apply(P"ZZI", C4) == P"+ZZZ"
  @test apply(P"IYZ", C4) == P"+YYX"

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
  C = one(CliffordOperator, 3)
  reduced, phase = reducetoX(P"XYZ", 2)
  @test reduced == P"XXZ"
  apply!(C, phase)
  @test C == C"XII IYI IIX ZII IZI IIZ" # apply! for Cliffords ignores phase

  C = one(CliffordOperator, 3)
  reduced, phase = reducetoX(P"YYY", 1)
  @test reduced == P"XYY"
  apply!(C, phase)
  @test C == C"YII IXI IIX ZII IZI IIZ"
end

@testset "build_D" begin
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

  test_D = one(CliffordOperator, 5)
  apply!(test_D, tCY, [3,1])
  apply!(test_D, tCX, [3,2])
  apply!(test_D, tCZ, [3,4])
  apply!(test_D, tCX, [3,5])
  @test build_D(P"YXXZX", 3) == test_D
end

@testset "disentangler" begin
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

  Dtest = one(CliffordOperator, 3)
  apply!(Dtest, tCY, [2,1])
  C = C"IXI XII IIX IZI ZII IIZ"
  P = P"XYI"
  @assert paulinature(1, C, P) == :disentanglable
  @test disentangler(1, C, P)[1] == Dtest

  Dtest = one(CliffordOperator, 3)
  apply!(Dtest, tPhase, [2]) # QuantumClifford's apply!(CliffordOperator, CliffordOperator) ignores phases
  apply!(Dtest, tCY, [2,1])
  C = C"IXI XII IIX IZI ZII IIZ"
  P = P"YYI"
  @assert paulinature(1, C, P) == :disentanglable
  @test disentangler(1, C, P)[1] == Dtest

  # No swap in disentangler: build_D keyed at i = 3 (the original free X position).
  Dtest = one(CliffordOperator, 3)
  apply!(Dtest, tCY, [3,1])
  C = C"IXI XII IIX IZI ZII IIZ"
  P = P"IYX"
  @assert paulinature(1, C, P) == :disentanglable
  @test disentangler(1, C, P)[1] == Dtest
end
