# Golden analytic cases for `addmagicstate!` and `apply!` (T3/T5). See test/CONVENTIONS.md.
# Every expected value is derived by hand in the comments, not from the code.
#
# Conventions used (CONVENTIONS.md §2, §4):
#   gate (P, ϕ) = exp(iϕP) = cos ϕ + i sin ϕ · P, for every signed ϕ;
#   old frame: |ψ⟩ = Cdag⁻¹ |mps⟩;
#   addmagicstate!(ψ, k, φ) multiplies qubit k+1 of the MPS by cos φ + i sin φ · X.
#
# Single-qubit facts. For a|0⟩ + b|1⟩: ⟨X⟩ = 2 Re(a* b), ⟨Y⟩ = 2 Im(a* b), ⟨Z⟩ = |a|² − |b|².
#   exp(iφX)|0⟩ = cos φ|0⟩ + i sin φ|1⟩   →  (⟨X⟩, ⟨Y⟩, ⟨Z⟩) = (0, sin 2φ, cos 2φ)
#   exp(iφY)|0⟩ = cos φ|0⟩ −   sin φ|1⟩   →  (−sin 2φ, 0, cos 2φ)
#   exp(iφZ)|0⟩ = e^{iφ}|0⟩               →  (0, 0, 1)

using DisentangleCAMPS: addmagicstate!
using CliffordMPS
using QuantumClifford
using Test

"Real expectation value; the imaginary part must vanish."
function expval(ψ, P)
  ev = expectation(ψ, P)
  @test abs(imag(ev)) <= 1e-10
  return real(ev)
end

camps_bloch(ψ) = (expval(ψ, P"X"), expval(ψ, P"Y"), expval(ψ, P"Z"))

function test_bloch(ψ, x, y, z)
  bx, by, bz = camps_bloch(ψ)
  @test bx ≈ x atol = 1e-10
  @test by ≈ y atol = 1e-10
  @test bz ≈ z atol = 1e-10
end

"`gate` on qubits `inds` of an N-qubit register, as a CliffordOperator."
function embed(gate, N, inds)
  D = one(CliffordOperator, N)
  apply!(D, gate, inds)
  return D
end

# Angles away from every Clifford angle, plus Clifford ones where the formula still holds.
const MAGIC_ANGLES = [0.0, π/8, π/4, 3π/8, π/3, π/2, -π/4, -π/8, 1.1, 2.3, -0.7]

@testset "addmagicstate!" begin
  @testset "return value" begin
    @test addmagicstate!(CAMPS(1), 0, 0.3) === nothing
  end

  @testset "N=1, identity frame" begin
    # cos φ|0⟩ + i sin φ|1⟩ → (0, sin 2φ, cos 2φ)
    for φ in MAGIC_ANGLES
      ψ = CAMPS(1)
      addmagicstate!(ψ, 0, φ)
      test_bloch(ψ, 0.0, sin(2φ), cos(2φ))
    end
  end

  @testset "acts on qubit k+1 only, N=3" begin
    φ = π/3   # (0, √3/2, −1/2) on the magic qubit; (0, 0, 1) on the others
    for k in 0:2
      ψ = CAMPS(3)
      addmagicstate!(ψ, k, φ)
      for j in 1:3
        x, y, z = j == k + 1 ? (0.0, sqrt(3)/2, -0.5) : (0.0, 0.0, 1.0)
        @test expval(ψ, single_x(3, j)) ≈ x atol = 1e-10
        @test expval(ψ, single_y(3, j)) ≈ y atol = 1e-10
        @test expval(ψ, single_z(3, j)) ≈ z atol = 1e-10
      end
    end
  end

  @testset "two magic qubits, product state, N=2" begin
    # Qubit 1 with φ₁ = π/8: (0, √2/2, √2/2). Qubit 2 with φ₂ = π/3: (0, √3/2, −1/2).
    # Product state: ⟨P₁P₂⟩ = ⟨P₁⟩⟨P₂⟩.
    ψ = CAMPS(2)
    addmagicstate!(ψ, 0, π/8)
    addmagicstate!(ψ, 1, π/3)
    @test expval(ψ, P"YY") ≈  sqrt(2)/2 * sqrt(3)/2 atol = 1e-10
    @test expval(ψ, P"ZZ") ≈ -sqrt(2)/4             atol = 1e-10
    @test expval(ψ, P"YZ") ≈ -sqrt(2)/4             atol = 1e-10
    @test expval(ψ, P"ZY") ≈  sqrt(2)/2 * sqrt(3)/2 atol = 1e-10
    @test expval(ψ, P"YI") ≈  sqrt(2)/2             atol = 1e-10
    @test expval(ψ, P"IZ") ≈ -0.5                   atol = 1e-10
    @test expval(ψ, P"XX") ≈  0.0                   atol = 1e-10
    @test expval(ψ, P"XY") ≈  0.0                   atol = 1e-10
  end

  # Non-identity frames (frame-dependent: these construct `Cdag` directly).
  @testset "Cdag = S, N=1" begin
    # |ψ⟩ = S†(cos φ|0⟩ + i sin φ|1⟩) = cos φ|0⟩ + sin φ|1⟩ → (sin 2φ, 0, cos 2φ)
    for φ in [π/8, π/4, π/3, -0.7]
      ψ = CAMPS(embed(tPhase, 1, [1]))
      addmagicstate!(ψ, 0, φ)
      test_bloch(ψ, sin(2φ), 0.0, cos(2φ))
    end
  end

  @testset "Cdag = H, N=1" begin
    # |ψ⟩ = H(cos φ|0⟩ + i sin φ|1⟩) = (e^{iφ}|0⟩ + e^{−iφ}|1⟩)/√2 → (cos 2φ, −sin 2φ, 0)
    for φ in [π/8, π/4, π/3, -0.7]
      ψ = CAMPS(embed(tHadamard, 1, [1]))
      addmagicstate!(ψ, 0, φ)
      test_bloch(ψ, cos(2φ), -sin(2φ), 0.0)
    end
  end

  @testset "Cdag = CNOT(1→2), N=2" begin
    # |ψ⟩ = CNOT (cos φ|0⟩ + i sin φ|1⟩)|0⟩ = cos φ|00⟩ + i sin φ|11⟩ =: a|00⟩ + b|11⟩.
    #   ⟨ZI⟩ = ⟨IZ⟩ = cos 2φ, ⟨ZZ⟩ = 1.
    #   XX and YY swap |00⟩ ↔ |11⟩ with real factors: ⟨XX⟩ = 2 Re(a* b) = 0, ⟨YY⟩ = −2 Re(a* b) = 0.
    #   XY|00⟩ = i|11⟩, XY|11⟩ = −i|00⟩: ⟨XY⟩ = a*(−i b) + b*(i a) = 2 Im(a* b) = sin 2φ. Same for YX.
    for φ in [π/8, π/4, π/3, -0.7]
      ψ = CAMPS(embed(tCNOT, 2, [1, 2]))
      addmagicstate!(ψ, 0, φ)
      @test expval(ψ, P"ZI") ≈ cos(2φ) atol = 1e-10
      @test expval(ψ, P"IZ") ≈ cos(2φ) atol = 1e-10
      @test expval(ψ, P"ZZ") ≈ 1.0     atol = 1e-10
      @test expval(ψ, P"XX") ≈ 0.0     atol = 1e-10
      @test expval(ψ, P"YY") ≈ 0.0     atol = 1e-10
      @test expval(ψ, P"XY") ≈ sin(2φ) atol = 1e-10
      @test expval(ψ, P"YX") ≈ sin(2φ) atol = 1e-10
      @test expval(ψ, P"XI") ≈ 0.0     atol = 1e-10
    end
  end
end

@testset "apply!(ψ, k, P, ϕ): golden cases" begin
  @testset "N=1, general angles" begin
    for ϕ in [π/8, π/3, 1.1, 2.3, -0.7, -π/8]
      ψ = CAMPS(1)
      k = apply!(ψ, 0, P"X", ϕ)
      @test k == 1
      test_bloch(ψ, 0.0, sin(2ϕ), cos(2ϕ))

      ψ = CAMPS(1)
      k = apply!(ψ, 0, P"Y", ϕ)
      @test k == 1
      test_bloch(ψ, -sin(2ϕ), 0.0, cos(2ϕ))

      ψ = CAMPS(1)
      k = apply!(ψ, 0, P"Z", ϕ)
      @test k == 0                       # only a phase on |0⟩: nothing to do
      test_bloch(ψ, 0.0, 0.0, 1.0)

      # Negative Pauli: exp(iϕ(−X)) = exp(−iϕX) → (0, −sin 2ϕ, cos 2ϕ)
      ψ = CAMPS(1)
      apply!(ψ, 0, P"-X", ϕ)
      test_bloch(ψ, 0.0, -sin(2ϕ), cos(2ϕ))
    end
  end

  @testset "N=1, two rotations about X add up" begin
    # exp(iϕ₂X) exp(iϕ₁X)|0⟩ = exp(i(ϕ₁+ϕ₂)X)|0⟩; the second gate is :logical (k = 1).
    ϕ1, ϕ2 = 0.4, 0.9
    ψ = CAMPS(1)
    k = apply!(ψ, 0, P"X", ϕ1)
    k = apply!(ψ, k, P"X", ϕ2)
    @test k == 1
    test_bloch(ψ, 0.0, sin(2 * (ϕ1 + ϕ2)), cos(2 * (ϕ1 + ϕ2)))
  end

  @testset "N=2, exp(iϕ XX)|00⟩" begin
    # cos ϕ|00⟩ + i sin ϕ|11⟩: the same state as in the CNOT-frame case above.
    for ϕ in [π/8, π/3, -0.7, π/4]         # π/4 takes the Clifford branch: (|00⟩ + i|11⟩)/√2
      ψ = CAMPS(2)
      apply!(ψ, 0, P"XX", ϕ)
      @test expval(ψ, P"ZI") ≈ cos(2ϕ) atol = 1e-10
      @test expval(ψ, P"IZ") ≈ cos(2ϕ) atol = 1e-10
      @test expval(ψ, P"ZZ") ≈ 1.0     atol = 1e-10
      @test expval(ψ, P"XX") ≈ 0.0     atol = 1e-10
      @test expval(ψ, P"YY") ≈ 0.0     atol = 1e-10
      @test expval(ψ, P"XY") ≈ sin(2ϕ) atol = 1e-10
      @test expval(ψ, P"YX") ≈ sin(2ϕ) atol = 1e-10
    end
  end

  @testset "N=1, Clifford angles" begin
    # exp(iϕX)|0⟩ → (0, sin 2ϕ, cos 2ϕ) holds for every ϕ:
    #   ϕ = 0, π:   (0, 0, 1)          ϕ = ±π/2: (0, 0, −1)   (±iX|0⟩ = ±i|1⟩)
    #   ϕ = +π/4:   (0, +1, 0)         ϕ = −π/4: (0, −1, 0)
    #   ϕ = +3π/4:  (0, −1, 0)         ϕ = −3π/4: (0, +1, 0)
    function xrot(ϕ)
      ψ = CAMPS(1)
      k = apply!(ψ, 0, P"X", ϕ)
      return k, camps_bloch(ψ)
    end
    close3(b, t) = all(isapprox.(b, t; atol = 1e-10))

    for (ϕ, target) in [(0.0, (0, 0, 1)), (π, (0, 0, 1)), (-π, (0, 0, 1)),
                        (π/2, (0, 0, -1)), (-π/2, (0, 0, -1)),
                        (π/4, (0, 1, 0)), (-3π/4, (0, 1, 0))]
      k, b = xrot(ϕ)
      @test k == 0                       # Clifford gates add no magic qubit
      @test close3(b, target)
    end

    # Known bugs in the old `apply!` (plan, hard spot 6): it branches on abs(ϕ).
    # These must turn green after Phase 2.
    @test_broken close3(xrot(-π/4)[2], (0, -1, 0))
    @test_broken close3(xrot(3π/4)[2], (0, -1, 0))
  end

  @testset "N=1, ±π/2 after another gate" begin
    close3(b, t) = all(isapprox.(b, t; atol = 1e-10))
    function two_gates(P1, ϕ1, P2, ϕ2)
      ψ = CAMPS(1)
      k = apply!(ψ, 0, P1, ϕ1)
      apply!(ψ, k, P2, ϕ2)
      return camps_bloch(ψ)
    end

    for ϕ in (π/2, -π/2)
      # General X rotation, then ±iZ: (0, sin 2ϕ₁, cos 2ϕ₁) → Z flips ⟨X⟩ and ⟨Y⟩.
      @test close3(two_gates(P"X", 0.4, P"Z", ϕ), (0.0, -sin(0.8), cos(0.8)))

      # Known bug in the old `apply!` (plan, hard spot 6): the ±π/2 branch multiplies the
      # frame on the wrong side. These must turn green after Phase 2.
      #
      # exp(iπ/4·X)|0⟩ = (|0⟩ + i|1⟩)/√2, then ±iZ → (|0⟩ − i|1⟩)/√2 up to a phase: (0, −1, 0).
      @test_broken close3(two_gates(P"X", π/4, P"Z", ϕ), (0.0, -1.0, 0.0))
      # General Y rotation, cos ϕ₁|0⟩ − sin ϕ₁|1⟩ = (−sin 2ϕ₁, 0, cos 2ϕ₁), then ±iX:
      # X keeps ⟨X⟩ and flips ⟨Y⟩, ⟨Z⟩ → (−sin 2ϕ₁, 0, −cos 2ϕ₁).
      @test_broken close3(two_gates(P"Y", 0.4, P"X", ϕ), (-sin(0.8), 0.0, -cos(0.8)))
    end
  end
end
