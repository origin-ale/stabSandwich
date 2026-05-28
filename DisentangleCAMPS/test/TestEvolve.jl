using DisentangleCAMPS: addmagicstate!
using CliffordMPS
using QuantumClifford
using Test

# Pauli expectation values are real; extract the real part from the ComplexF64 return.
expval(ψ, args...) = (expectation(ψ, args...))

@testset "addmagicstate! return value" begin
    ψ = CAMPS(1)
    @test addmagicstate!(ψ, 1, "X", 0.0) === nothing
    ψ = CAMPS(1)
    @test addmagicstate!(ψ, 1, "Y", π/4) === nothing
end

@testset "addmagicstate! op=X, N=1" begin
    # Operator applied to |0⟩: M = cos(θ)·I - i·sin(θ)·X
    # Resulting state: cos(θ)|0⟩ - i·sin(θ)|1⟩
    # Bloch vector: (⟨X⟩, ⟨Y⟩, ⟨Z⟩) = (0, -sin(2θ), cos(2θ))  →  Y-Z plane
    for (θ, sin2θ, cos2θ) in [
        (0.0,              0.0,          1.0          ),  # |0⟩
        (π/8,   sqrt(2)/2,    sqrt(2)/2  ),  # sin(π/4), cos(π/4)
        (π/4,              1.0,          0.0          ),  # sin(π/2)=1, cos(π/2)=0
        (3π/8,  sqrt(2)/2,   -sqrt(2)/2  ),  # sin(3π/4), cos(3π/4)
        (π/3,   sqrt(3)/2,   -0.5        ),  # sin(2π/3), cos(2π/3)
        (π/2,              0.0,         -1.0          ),  # |1⟩ up to global phase
        (-π/4,            -1.0,          0.0          ),  # sin(-π/2), cos(-π/2)
    ]
        ψ = CAMPS(1)
        addmagicstate!(ψ, 1, "X", θ)
        @test expval(ψ, P"X") ≈  0.0    atol=1e-10
        @test expval(ψ, P"Y") ≈ -sin2θ  atol=1e-10
        @test expval(ψ, P"Z") ≈  cos2θ  atol=1e-10
    end
end

@testset "addmagicstate! op=Y, N=1" begin
    # Operator applied to |0⟩: M = cos(θ)·I + sin(θ)·X
    # Resulting state: cos(θ)|0⟩ + sin(θ)|1⟩
    # Bloch vector: (⟨X⟩, ⟨Y⟩, ⟨Z⟩) = (sin(2θ), 0, cos(2θ))  →  X-Z plane
    for (θ, sin2θ, cos2θ) in [
        (0.0,              0.0,          1.0          ),  # |0⟩
        (π/8,   sqrt(2)/2,    sqrt(2)/2  ),  # sin(π/4), cos(π/4)
        (π/4,              1.0,          0.0          ),  # |+⟩
        (3π/8,  sqrt(2)/2,   -sqrt(2)/2  ),
        (π/3,   sqrt(3)/2,   -0.5        ),
        (π/2,              0.0,         -1.0          ),  # |1⟩
        (-π/4,            -1.0,          0.0          ),  # |−⟩
    ]
        ψ = CAMPS(1)
        addmagicstate!(ψ, 1, "Y", θ)
        @test expval(ψ, P"X") ≈  sin2θ  atol=1e-10
        @test expval(ψ, P"Y") ≈  0.0    atol=1e-10
        @test expval(ψ, P"Z") ≈  cos2θ  atol=1e-10
    end
end

@testset "addmagicstate! Bloch sphere purity" begin
    # Any pure single-qubit state satisfies ⟨X⟩² + ⟨Y⟩² + ⟨Z⟩² = 1.
    # Tested at representative angles, including non-special values.
    for θ in [0.0, π/8, π/4, π/3, π/2, -π/4, 1.1, 2.3]
        for op in ("X", "Y")
            ψ = CAMPS(1)
            addmagicstate!(ψ, 1, op, θ)
            x, y, z = expval(ψ, P"X"), expval(ψ, P"Y"), expval(ψ, P"Z")
            @test x^2 + y^2 + z^2 ≈ 1.0 atol=1e-10
        end
    end
end

@testset "addmagicstate! single-site expectations in N=3" begin
    # The magic qubit carries the expected Bloch-vector components.
    # Untouched qubits remain in |0⟩: ⟨X⟩=⟨Y⟩=0, ⟨Z⟩=1.
    for i in 1:3
        for (op, θ) in [("X", π/4), ("Y", π/4), ("X", π/8), ("Y", π/3)]
            ψ = CAMPS(3)
            addmagicstate!(ψ, i, op, θ)
            sin2θ, cos2θ = sin(2θ), cos(2θ)

            if op == "X"
                @test expval(ψ, P"X", [i]) ≈  0.0    atol=1e-10
                @test expval(ψ, P"Y", [i]) ≈ -sin2θ  atol=1e-10
                @test expval(ψ, P"Z", [i]) ≈  cos2θ  atol=1e-10
            else
                @test expval(ψ, P"X", [i]) ≈  sin2θ  atol=1e-10
                @test expval(ψ, P"Y", [i]) ≈  0.0    atol=1e-10
                @test expval(ψ, P"Z", [i]) ≈  cos2θ  atol=1e-10
            end

            for j in setdiff(1:3, [i])
                @test expval(ψ, P"X", [j]) ≈ 0.0 atol=1e-10
                @test expval(ψ, P"Y", [j]) ≈ 0.0 atol=1e-10
                @test expval(ψ, P"Z", [j]) ≈ 1.0 atol=1e-10
            end
        end
    end
end

@testset "addmagicstate! multi-qubit Pauli strings in N=3" begin
    # For a product state ⟨P₁ P₂ P₃⟩ = ⟨P₁⟩⟨P₂⟩⟨P₃⟩.
    # Untouched qubits in |0⟩ contribute: ⟨Z⟩=1, ⟨X⟩=⟨Y⟩=0, ⟨I⟩=1.

    # Qubit 2: op="Y", θ=π/4  →  |+⟩,  ⟨X₂⟩=1, ⟨Y₂⟩=0, ⟨Z₂⟩=0
    ψ = CAMPS(3); addmagicstate!(ψ, 2, "Y", π/4)
    @test expval(ψ, P"ZXZ") ≈  1.0 atol=1e-10   # 1·1·1
    @test expval(ψ, P"ZYZ") ≈  0.0 atol=1e-10   # 1·0·1
    @test expval(ψ, P"ZZZ") ≈  0.0 atol=1e-10   # 1·0·1   (⟨Z₂⟩=0 at θ=π/4)
    @test expval(ψ, P"XXI") ≈  0.0 atol=1e-10   # 0·1·1   (⟨X₁⟩=0 on |0⟩)
    @test expval(ψ, P"IXI") ≈  1.0 atol=1e-10   # 1·1·1

    # Qubit 2: op="X", θ=π/4  →  (|0⟩-i|1⟩)/√2,  ⟨X₂⟩=0, ⟨Y₂⟩=-1, ⟨Z₂⟩=0
    ψ = CAMPS(3); addmagicstate!(ψ, 2, "X", π/4)
    @test expval(ψ, P"ZYZ") ≈ -1.0 atol=1e-10   # 1·(-1)·1
    @test expval(ψ, P"IXI") ≈  0.0 atol=1e-10   # 1·0·1
    @test expval(ψ, P"IZI") ≈  0.0 atol=1e-10   # 1·0·1   (⟨Z₂⟩=0)
    @test expval(ψ, P"ZZZ") ≈  0.0 atol=1e-10   # 1·0·1

    # Qubit 1: op="X", θ=π/8  →  ⟨Y₁⟩=-√2/2, ⟨Z₁⟩=√2/2, ⟨X₁⟩=0
    ψ = CAMPS(3); addmagicstate!(ψ, 1, "X", π/8)
    @test expval(ψ, P"YZZ") ≈ -sqrt(2)/2 atol=1e-10  # (-√2/2)·1·1
    @test expval(ψ, P"ZZZ") ≈  sqrt(2)/2 atol=1e-10  # (√2/2)·1·1
    @test expval(ψ, P"XII") ≈  0.0       atol=1e-10  # 0·1·1   (⟨X₁⟩=0 for op="X")

    # Qubit 3: op="Y", θ=π/3  →  ⟨X₃⟩=√3/2, ⟨Y₃⟩=0, ⟨Z₃⟩=-1/2
    ψ = CAMPS(3); addmagicstate!(ψ, 3, "Y", π/3)
    @test expval(ψ, P"ZZX") ≈  sqrt(3)/2 atol=1e-10  # 1·1·(√3/2)
    @test expval(ψ, P"ZZZ") ≈ -0.5       atol=1e-10  # 1·1·(-1/2)
    @test expval(ψ, P"IIY") ≈  0.0       atol=1e-10  # 1·1·0   (⟨Y₃⟩=0 for op="Y")
end

@testset "addmagicstate! two qubits made magic, N=2" begin
    # Qubit 1: op="X", θ=π/4  →  ⟨X₁⟩=0,    ⟨Y₁⟩=-1,    ⟨Z₁⟩=0
    # Qubit 2: op="Y", θ=π/4  →  ⟨X₂⟩=1,    ⟨Y₂⟩=0,     ⟨Z₂⟩=0
    ψ = CAMPS(2)
    addmagicstate!(ψ, 1, "X", π/4)
    addmagicstate!(ψ, 2, "Y", π/4)
    @test expval(ψ, P"YX") ≈ -1.0 atol=1e-10  # (-1)·1
    @test expval(ψ, P"ZZ") ≈  0.0 atol=1e-10  # 0·0
    @test expval(ψ, P"ZX") ≈  0.0 atol=1e-10  # 0·1
    @test expval(ψ, P"YZ") ≈  0.0 atol=1e-10  # (-1)·0
    @test expval(ψ, P"IX") ≈  1.0 atol=1e-10  # 1·1
    @test expval(ψ, P"YI") ≈ -1.0 atol=1e-10  # (-1)·1

    # Qubit 1: op="Y", θ=π/8  →  ⟨X₁⟩=√2/2,  ⟨Y₁⟩=0,     ⟨Z₁⟩=√2/2
    # Qubit 2: op="X", θ=π/8  →  ⟨X₂⟩=0,    ⟨Y₂⟩=-√2/2, ⟨Z₂⟩=√2/2
    ψ = CAMPS(2)
    addmagicstate!(ψ, 1, "Y", π/8)
    addmagicstate!(ψ, 2, "X", π/8)
    @test expval(ψ, P"XY") ≈ -0.5 atol=1e-10  # (√2/2)·(-√2/2) = -1/2
    @test expval(ψ, P"ZZ") ≈  0.5 atol=1e-10  # (√2/2)·(√2/2) = 1/2
    @test expval(ψ, P"XZ") ≈  0.5 atol=1e-10  # (√2/2)·(√2/2) = 1/2
    @test expval(ψ, P"ZY") ≈ -0.5 atol=1e-10  # (√2/2)·(-√2/2) = -1/2
    @test expval(ψ, P"XI") ≈  sqrt(2)/2 atol=1e-10  # √2/2·1
    @test expval(ψ, P"IY") ≈ -sqrt(2)/2 atol=1e-10  # 1·(-√2/2)
end

# ── Non-identity Clifford frame ────────────────────────────────────────────────
#
# The expectation value of P on a CAMPS(Cdag) state after addmagicstate!(ψ,i,op,θ)
# is ⟨mps|(P^Cdag)|mps⟩ = ⟨mps|CdagPCdag†|mps⟩, where the physical state is
# inv(Cdag)|mps⟩.  Conjugating by Cdag rotates the Pauli basis, so the
# numerically observed single-qubit expectations differ from the identity case.
# All literals below are derived independently by applying the Clifford action
# on the physical qubit state.

@testset "addmagicstate! Cdag=Phase(S), N=1" begin
    # S†|0⟩=|0⟩, S†|1⟩=−i|1⟩.
    # op="X": MPS cos(θ)|0⟩−i sin(θ)|1⟩  →  physical cos(θ)|0⟩−sin(θ)|1⟩
    #   ⟨X⟩=−sin(2θ),  ⟨Y⟩=0,      ⟨Z⟩=cos(2θ)
    # op="Y": MPS cos(θ)|0⟩+sin(θ)|1⟩    →  physical cos(θ)|0⟩−i sin(θ)|1⟩
    #   ⟨X⟩=0,         ⟨Y⟩=−sin(2θ), ⟨Z⟩=cos(2θ)
    Cs = one(CliffordOperator, 1)
    apply!(Cs, tPhase, [1])

    for (θ, sin2θ, cos2θ) in [
        (π/8, sqrt(2)/2,  sqrt(2)/2),
        (π/4, 1.0,        0.0      ),
        (π/3, sqrt(3)/2, -0.5      ),
    ]
        ψ = CAMPS(Cs); addmagicstate!(ψ, 1, "X", θ)
        @test expval(ψ, P"X") ≈ -sin2θ  atol=1e-10
        @test expval(ψ, P"Y") ≈  0.0    atol=1e-10
        @test expval(ψ, P"Z") ≈  cos2θ  atol=1e-10

        ψ = CAMPS(Cs); addmagicstate!(ψ, 1, "Y", θ)
        @test expval(ψ, P"X") ≈  0.0    atol=1e-10
        @test expval(ψ, P"Y") ≈ -sin2θ  atol=1e-10
        @test expval(ψ, P"Z") ≈  cos2θ  atol=1e-10
    end
end

@testset "addmagicstate! Cdag=Hadamard, N=1" begin
    # H is self-inverse, so physical = H|mps⟩.
    # op="X": H(cos θ|0⟩−i sin θ|1⟩) = [e^{−iθ}|0⟩+e^{iθ}|1⟩]/√2
    #   ⟨X⟩=cos(2θ),  ⟨Y⟩=sin(2θ),  ⟨Z⟩=0
    # op="Y": H(cos θ|0⟩+sin θ|1⟩) = [(cos θ+sin θ)|0⟩+(cos θ−sin θ)|1⟩]/√2
    #   ⟨X⟩=cos(2θ),  ⟨Y⟩=0,        ⟨Z⟩=sin(2θ)
    Ch = one(CliffordOperator, 1)
    apply!(Ch, tHadamard, [1])

    for (θ, sin2θ, cos2θ) in [
        (π/8, sqrt(2)/2,  sqrt(2)/2),
        (π/4, 1.0,        0.0      ),
        (π/3, sqrt(3)/2, -0.5      ),
    ]
        ψ = CAMPS(Ch); addmagicstate!(ψ, 1, "X", θ)
        @test expval(ψ, P"X") ≈  cos2θ  atol=1e-10
        @test expval(ψ, P"Y") ≈  sin2θ  atol=1e-10
        @test expval(ψ, P"Z") ≈  0.0    atol=1e-10

        ψ = CAMPS(Ch); addmagicstate!(ψ, 1, "Y", θ)
        @test expval(ψ, P"X") ≈  cos2θ  atol=1e-10
        @test expval(ψ, P"Y") ≈  0.0    atol=1e-10
        @test expval(ψ, P"Z") ≈  sin2θ  atol=1e-10
    end
end

@testset "addmagicstate! Cdag=CNOT(1→2), N=2" begin
    # CNOT is self-inverse, so physical = CNOT|mps⟩.
    Cc = one(CliffordOperator, 2)
    apply!(Cc, tCNOT, [1, 2])

    # qubit 1, op="Y", θ=π/4: MPS = |+⟩⊗|0⟩
    # physical = CNOT(|+⟩⊗|0⟩) = (|00⟩+|11⟩)/√2  — Bell state |Φ+⟩
    # |Φ+⟩ correlators: ⟨ZZ⟩=1, ⟨XX⟩=1, ⟨YY⟩=−1; single-site: ⟨ZI⟩=⟨IZ⟩=⟨XI⟩=⟨IX⟩=0
    ψ = CAMPS(Cc); addmagicstate!(ψ, 1, "Y", π/4)
    @test expval(ψ, P"ZI") ≈  0.0 atol=1e-10
    @test expval(ψ, P"IZ") ≈  0.0 atol=1e-10
    @test expval(ψ, P"ZZ") ≈  1.0 atol=1e-10
    @test expval(ψ, P"XX") ≈  1.0 atol=1e-10
    @test expval(ψ, P"YY") ≈ -1.0 atol=1e-10
    @test expval(ψ, P"XI") ≈  0.0 atol=1e-10
    @test expval(ψ, P"IX") ≈  0.0 atol=1e-10

    # qubit 1, op="X", θ=π/4: MPS = (|0⟩−i|1⟩)/√2⊗|0⟩
    # physical = CNOT·MPS = (|00⟩−i|11⟩)/√2
    # ⟨ZZ⟩=1, ⟨XX⟩=0, ⟨YY⟩=0, ⟨ZI⟩=0
    ψ = CAMPS(Cc); addmagicstate!(ψ, 1, "X", π/4)
    @test expval(ψ, P"ZZ") ≈  1.0 atol=1e-10
    @test expval(ψ, P"XX") ≈  0.0 atol=1e-10
    @test expval(ψ, P"YY") ≈  0.0 atol=1e-10
    @test expval(ψ, P"ZI") ≈  0.0 atol=1e-10
    @test expval(ψ, P"IZ") ≈  0.0 atol=1e-10

    # qubit 1, op="X", θ=π/8: MPS = (cos(π/8)|0⟩−i sin(π/8)|1⟩)⊗|0⟩
    # physical = cos(π/8)|00⟩−i sin(π/8)|11⟩
    # ⟨ZI⟩ = cos(2·π/8) = cos(π/4) = √2/2  (Z on qubit 1 of physical state)
    ψ = CAMPS(Cc); addmagicstate!(ψ, 1, "X", π/8)
    @test expval(ψ, P"ZI") ≈  sqrt(2)/2 atol=1e-10
    @test expval(ψ, P"ZZ") ≈  1.0       atol=1e-10

    # qubit 2, op="Y", θ=π/4: MPS = |0⟩⊗|+⟩
    # physical = CNOT(|0⟩⊗|+⟩) = |0⟩⊗|+⟩  (control=0 → target unchanged)
    # ⟨ZI⟩=1, ⟨IX⟩=1, ⟨ZX⟩=1, ⟨ZZ⟩=0, ⟨XX⟩=0
    ψ = CAMPS(Cc); addmagicstate!(ψ, 2, "Y", π/4)
    @test expval(ψ, P"ZI") ≈  1.0 atol=1e-10
    @test expval(ψ, P"IX") ≈  1.0 atol=1e-10
    @test expval(ψ, P"ZX") ≈  1.0 atol=1e-10
    @test expval(ψ, P"ZZ") ≈  0.0 atol=1e-10
    @test expval(ψ, P"XX") ≈  0.0 atol=1e-10

    # qubit 2, op="X", θ=π/4: MPS = |0⟩⊗(|0⟩−i|1⟩)/√2
    # physical = |0⟩⊗(|0⟩−i|1⟩)/√2  (same — control=0)
    # ⟨ZI⟩=1, ⟨IY⟩=−1, ⟨ZY⟩=−1
    ψ = CAMPS(Cc); addmagicstate!(ψ, 2, "X", π/4)
    @test expval(ψ, P"ZI") ≈  1.0 atol=1e-10
    @test expval(ψ, P"IY") ≈ -1.0 atol=1e-10
    @test expval(ψ, P"ZY") ≈ -1.0 atol=1e-10

    # both qubits magical: qubit 1 op="Y" θ=π/4 → |+⟩, qubit 2 op="X" θ=π/4 → |−y⟩
    # MPS = |+⟩⊗|−y⟩  (⟨X₁⟩=1, ⟨Y₂⟩=−1, all others = 0)
    # physical = CNOT(|+⟩⊗|−y⟩) = (|00⟩−i|01⟩−i|10⟩+|11⟩)/2  — entangled
    # CNOT conjugation: XX→XI, ZY→IY, YZ→XY; all single-site expectations vanish
    ψ = CAMPS(Cc)
    addmagicstate!(ψ, 1, "Y", π/4)
    addmagicstate!(ψ, 2, "X", π/4)
    @test expval(ψ, P"ZI") ≈  0.0 atol=1e-10
    @test expval(ψ, P"IZ") ≈  0.0 atol=1e-10
    @test expval(ψ, P"XI") ≈  0.0 atol=1e-10
    @test expval(ψ, P"IX") ≈  0.0 atol=1e-10
    @test expval(ψ, P"XX") ≈  1.0 atol=1e-10
    @test expval(ψ, P"ZY") ≈ -1.0 atol=1e-10
    @test expval(ψ, P"YZ") ≈ -1.0 atol=1e-10
    @test expval(ψ, P"ZZ") ≈  0.0 atol=1e-10

    # both qubits magical: qubit 1 op="Y" θ=π/8, qubit 2 op="Y" θ=π/8
    # each MPS qubit: cos(π/8)|0⟩+sin(π/8)|1⟩  →  ⟨X⟩=√2/2, ⟨Y⟩=0, ⟨Z⟩=√2/2
    # physical = CNOT|m⟩⊗|m⟩  — entangled
    # CNOT conjugation: ZI→ZI, IZ→ZZ, ZZ→IZ, XX→XI, XI→XX, IX→IX, ZX→ZX, YY→−XZ
    ψ = CAMPS(Cc)
    addmagicstate!(ψ, 1, "Y", π/8)
    addmagicstate!(ψ, 2, "Y", π/8)
    @test expval(ψ, P"ZI") ≈  sqrt(2)/2 atol=1e-10
    @test expval(ψ, P"IZ") ≈  0.5       atol=1e-10
    @test expval(ψ, P"ZZ") ≈  sqrt(2)/2 atol=1e-10
    @test expval(ψ, P"XX") ≈  sqrt(2)/2 atol=1e-10
    @test expval(ψ, P"XI") ≈  0.5       atol=1e-10
    @test expval(ψ, P"IX") ≈  sqrt(2)/2 atol=1e-10
    @test expval(ψ, P"ZX") ≈  0.5       atol=1e-10
    @test expval(ψ, P"YY") ≈ -0.5       atol=1e-10
end

@testset "addmagicstate! Cdag=CY(1→2)·CX(2→1), |m⟩|m⟩" begin
    tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
    Ccyx = one(CliffordOperator, 2)
    apply!(Ccyx, tCY,   [1, 2])
    apply!(Ccyx, tCNOT, [2, 1])

    ψ = CAMPS(Ccyx)
    addmagicstate!(ψ, 1, "Y", +1, π/8)
    addmagicstate!(ψ, 2, "X", -1, π/8)

    @test expval(ψ, P"ZI") ≈  0.5         atol=1e-10
    @test expval(ψ, P"ZZ") ≈  sqrt(2)/2   atol=1e-10
    @test expval(ψ, P"XX") ≈  0           atol=1e-10
    @test expval(ψ, P"YY") ≈  0.0         atol=1e-10
    @test expval(ψ, P"IY") ≈  0.5         atol=1e-10
end
