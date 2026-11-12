# Golden set (T3): hand-verified expected values. See test/CONVENTIONS.md.
# Nothing here is derived from the implementation under test.

using QuantumClifford
using CliffordMPS   # defines `CliffordOperator * PauliOperator`
using Test

# Assumption pins: these test QuantumClifford (and one thin CliffordMPS wrapper),
# not our code. They fix the conventions of CONVENTIONS.md §4, which every other
# expected value relies on.
#
# Hand derivation, with S = diag(1, i) and H the Hadamard:
#   S X S† = +Y    S Y S† = -X    S Z S† = Z      (S† X S = -Y,  S† Y S = +X)
#   H X H  = +Z    H Y H  = -Y    H Z H  = X
#   CNOT(1→2): XI → XX, IX → IX, ZI → ZI, IZ → ZZ
# A CliffordOperator is defined by its tableau, the images C X_j C†, C Z_j C†.
@testset "assumption pins: Clifford conjugation direction" begin
  # The named gates carry the tableaus derived above.
  @test tPhase == C"Y Z"
  @test tHadamard == C"Z X"
  @test tCNOT == C"XX IX ZI ZZ"

  # C"…" lists the images of X_1…X_N first, then Z_1…Z_N.
  C = C"Y Z"
  @test C * Stabilizer([P"X"]) == Stabilizer([P"Y"])
  @test C * Stabilizer([P"Z"]) == Stabilizer([P"Z"])

  @testset "C * Stabilizer gives C P C†" begin
    @test (tPhase * Stabilizer([P"X"]))[1] == P"+Y"
    @test (tPhase * Stabilizer([P"Y"]))[1] == P"-X"
    @test (tPhase * Stabilizer([P"Z"]))[1] == P"+Z"
  end

  @testset "C * PauliOperator gives C P C†" begin
    @test tPhase * P"X" == P"+Y"
    @test tPhase * P"Y" == P"-X"
    @test tPhase * P"Z" == P"+Z"
    @test tHadamard * P"Y" == P"-Y"
    @test tCNOT * P"XI" == P"+XX"
    @test tCNOT * P"IZ" == P"+ZZ"
  end

  @testset "apply!(stabilizer, C) gives C P C†" begin
    @test apply!(Stabilizer([P"X"]), tPhase)[1] == P"+Y"
    @test apply!(Stabilizer([P"Y"]), tPhase)[1] == P"-X"
    @test apply!(Stabilizer([P"Y"]), tHadamard)[1] == P"-Y"
    @test apply!(Stabilizer([P"XI"]), tCNOT)[1] == P"+XX"
    @test apply!(Stabilizer([P"IZ"]), tCNOT)[1] == P"+ZZ"
  end

  @testset "inv(C) is C†: inv(C) * P and apply!(stab, inv(C)) give C† P C" begin
    @test inv(tPhase) * P"X" == P"-Y"
    @test inv(tPhase) * P"Y" == P"+X"
    @test apply!(Stabilizer([P"X"]), inv(tPhase))[1] == P"-Y"
  end

  @testset "A * B applies B first" begin
    # (HS) X (HS)† = H (S X S†) H = H Y H = -Y;  (SH) X (SH)† = S Z S† = Z.
    @test (tHadamard * tPhase) * P"X" == P"-Y"
    @test (tPhase * tHadamard) * P"X" == P"+Z"
    # (HS) Z (HS)† = H Z H = X;  (SH) Z (SH)† = S X S† = Y.
    @test (tHadamard * tPhase) * P"Z" == P"+X"
    @test (tPhase * tHadamard) * P"Z" == P"+Y"
  end
end
