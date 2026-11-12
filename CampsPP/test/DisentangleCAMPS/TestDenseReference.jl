# Tests of the dense reference simulator against hand-derived facts only.
# Run with: julia CampsPP/test/DisentangleCAMPS/TestDenseReference.jl
include("DenseReference.jl")
using .DenseReference
using Test

const ATOL = 1e-12
close(a, b) = all(isapprox.(a, b; atol = ATOL))
norm2(ψ) = sum(abs2, ψ)
# Matrix-vector product written out, so that the test needs nothing beyond Base.
matvec(A, v) = [sum(A[i, j] * v[j] for j in eachindex(v)) for i in axes(A, 1)]

# Literal single-qubit matrices in the basis (|0⟩, |1⟩).
const I2 = ComplexF64[1 0; 0 1]
const X2 = ComplexF64[0 1; 1 0]
const Y2 = ComplexF64[0 -im; im 0]
const Z2 = ComplexF64[1 0; 0 -1]

const ANGLES = [0.0, π/4, -π/4, π/2, -π/2, 3π/4, -3π/4, π, 0.3, -0.3, 1.1, -2.5, 7.0]

@testset "Dense reference" begin

    @testset "basis states and Pauli action" begin
        @test zero_state(1) == ComplexF64[1, 0]
        @test basis_state(1, [1]) == ComplexF64[0, 1]
        @test length(zero_state(5)) == 32
        @test zero_state(3) == basis_state(3, Int[])
        # X|0⟩ = |1⟩, X|1⟩ = |0⟩
        @test apply_pauli(ComplexF64[1, 0], "X") == ComplexF64[0, 1]
        @test apply_pauli(ComplexF64[0, 1], "X") == ComplexF64[1, 0]
        # Z|0⟩ = |0⟩, Z|1⟩ = -|1⟩
        @test apply_pauli(ComplexF64[1, 0], "Z") == ComplexF64[1, 0]
        @test apply_pauli(ComplexF64[0, 1], "Z") == ComplexF64[0, -1]
        # Y = [[0, -i], [i, 0]]: Y|0⟩ = i|1⟩ (first column), Y|1⟩ = -i|0⟩ (second column)
        @test apply_pauli(ComplexF64[1, 0], "Y") == ComplexF64[0, im]
        @test apply_pauli(ComplexF64[0, 1], "Y") == ComplexF64[-im, 0]
        # signs
        @test apply_pauli(ComplexF64[1, 0], "-X") == ComplexF64[0, -1]
        @test apply_pauli(ComplexF64[1, 0], "+X") == ComplexF64[0, 1]
        # every single-qubit Pauli agrees with its literal matrix on a generic vector
        v = ComplexF64[0.6 + 0.1im, -0.3 + 0.7im]
        for (s, M) in (("I", I2), ("X", X2), ("Y", Y2), ("Z", Z2))
            @test close(apply_pauli(v, s), matvec(M, v))
        end
    end

    @testset "single-qubit rotations of |0⟩" begin
        for ϕ in ANGLES
            c, s = cos(ϕ), sin(ϕ)
            # exp(iϕX)|0⟩ = cos ϕ|0⟩ + i sin ϕ X|0⟩ = cos ϕ|0⟩ + i sin ϕ|1⟩.
            # ⟨X⟩ = 2 Re(conj(c)·is) = 0, ⟨Y⟩ = 2 Im(conj(c)·is) = 2cs = sin 2ϕ,
            # ⟨Z⟩ = c² − s² = cos 2ϕ.   (For a|0⟩+b|1⟩: ⟨X⟩ = 2Re(āb), ⟨Y⟩ = 2Im(āb).)
            ψ = rotate(zero_state(1), "X", ϕ)
            @test close(ψ, [c, im * s])
            @test close(collect(bloch(ψ, 1)), [0, sin(2ϕ), cos(2ϕ)])
            # exp(iϕY)|0⟩ = cos ϕ|0⟩ + i sin ϕ·(i|1⟩) = cos ϕ|0⟩ − sin ϕ|1⟩.
            # ⟨X⟩ = 2Re(c·(−s)) = −sin 2ϕ, ⟨Y⟩ = 0 (real amplitudes), ⟨Z⟩ = cos 2ϕ.
            ψ = rotate(zero_state(1), "Y", ϕ)
            @test close(ψ, [c, -s])
            @test close(collect(bloch(ψ, 1)), [-sin(2ϕ), 0, cos(2ϕ)])
            # exp(iϕZ)|0⟩ = (cos ϕ + i sin ϕ)|0⟩ = e^{iϕ}|0⟩: a global phase.
            ψ = rotate(zero_state(1), "Z", ϕ)
            @test close(ψ, [cis(ϕ), 0])
            @test close(collect(bloch(ψ, 1)), [0, 0, 1])
            # exp(iϕZ)|1⟩ = e^{−iϕ}|1⟩, Bloch (0, 0, −1).
            ψ = rotate(basis_state(1, [1]), "Z", ϕ)
            @test close(ψ, [0, cis(-ϕ)])
            @test close(collect(bloch(ψ, 1)), [0, 0, -1])
        end
        # Special values from CONVENTIONS §2: exp(±iπ/4 P) = (1 ± iP)/√2,
        # exp(±iπ/2 P) = ±iP, exp(i 3π/4 P) = (−1 + iP)/√2, ϕ = 0, π give ±1.
        v = ComplexF64[0.6 + 0.1im, -0.3 + 0.7im]
        for P in ("X", "Y", "Z", "-Y")
            Pv = apply_pauli(v, P)
            @test close(rotate(v, P, π/4), (v + im * Pv) / sqrt(2))
            @test close(rotate(v, P, -π/4), (v - im * Pv) / sqrt(2))
            @test close(rotate(v, P, π/2), im * Pv)
            @test close(rotate(v, P, -π/2), -im * Pv)
            @test close(rotate(v, P, 3π/4), (-v + im * Pv) / sqrt(2))
            @test close(rotate(v, P, 0), v)
            @test close(rotate(v, P, π), -v)
        end
    end

    @testset "sign of the Pauli string" begin
        # exp(iϕ(−X)) = cos ϕ − i sin ϕ X = exp(−iϕX); on |0⟩: cos ϕ|0⟩ − i sin ϕ|1⟩.
        v = ComplexF64[0.6 + 0.1im, -0.3 + 0.7im]
        for ϕ in ANGLES
            @test close(rotate(zero_state(1), "-X", ϕ), [cos(ϕ), -im * sin(ϕ)])
            @test close(rotate(v, "-X", ϕ), rotate(v, "X", -ϕ))
            @test close(rotate(v, "+Y", ϕ), rotate(v, "Y", ϕ))
        end
        @test close(apply_pauli(v, "-Z"), -apply_pauli(v, "Z"))
        @test isapprox(expect(v, "-Z"), -expect(v, "Z"); atol = ATOL)
    end

    @testset "rotate vs rotate!" begin
        ψ = zero_state(2)
        φ = rotate(ψ, "XY", 0.4)
        @test ψ == zero_state(2)            # rotate does not mutate
        @test rotate!(ψ, "XY", 0.4) === ψ   # rotate! does, and returns its argument
        @test ψ == φ
    end

    @testset "two qubits: exp(iϕ XX)|00⟩" begin
        for ϕ in ANGLES
            c, s = cos(ϕ), sin(ϕ)
            # XX|00⟩ = |11⟩, so ψ = c|00⟩ + i s|11⟩.
            ψ = rotate(zero_state(2), "XX", ϕ)
            @test close(ψ, c * basis_state(2, Int[]) + im * s * basis_state(2, [1, 2]))
            # ZI, IZ: +1 on |00⟩, −1 on |11⟩ → c² − s² = cos 2ϕ. ZZ: +1 on both → 1.
            @test isapprox(expect(ψ, "ZI"), cos(2ϕ); atol = ATOL)
            @test isapprox(expect(ψ, "IZ"), cos(2ϕ); atol = ATOL)
            @test isapprox(expect(ψ, "ZZ"), 1; atol = ATOL)
            # XXψ = c|11⟩ + i s|00⟩; ⟨ψ|XXψ⟩ = c·(is) + conj(is)·c = 0.
            @test isapprox(expect(ψ, "XX"), 0; atol = ATOL)
            # YY|00⟩ = (i)(i)|11⟩ = −|11⟩, YY|11⟩ = (−i)(−i)|00⟩ = −|00⟩;
            # YYψ = −c|11⟩ − i s|00⟩; ⟨ψ|YYψ⟩ = c·(−is) + conj(is)·(−c) = −ics + ics = 0.
            @test isapprox(expect(ψ, "YY"), 0; atol = ATOL)
            # XY|00⟩ = |1⟩ ⊗ i|1⟩ = i|11⟩, XY|11⟩ = |0⟩ ⊗ (−i)|0⟩ = −i|00⟩;
            # XYψ = ic|11⟩ + (is)(−i)|00⟩ = s|00⟩ + ic|11⟩;
            # ⟨ψ|XYψ⟩ = c·s + conj(is)·(ic) = cs + sc = sin 2ϕ. YX is the same by symmetry.
            @test isapprox(expect(ψ, "XY"), sin(2ϕ); atol = ATOL)
            @test isapprox(expect(ψ, "YX"), sin(2ϕ); atol = ATOL)
            # Single-qubit X, Y vanish (they connect |00⟩ only to |10⟩ or |01⟩).
            for P in ("XI", "IX", "YI", "IY")
                @test isapprox(expect(ψ, P), 0; atol = ATOL)
            end
        end
    end

    @testset "exp(iπ/4 (XX+YY+ZZ)) = e^{iπ/4} SWAP" begin
        # XX + YY + ZZ = 2·SWAP − 1 (it is +1 on the triplet, −3 on the singlet), and
        # SWAP² = 1, so exp(iπ/4(2 SWAP − 1)) = e^{−iπ/4}·exp(iπ/2 SWAP) = e^{−iπ/4}·i·SWAP
        # = e^{iπ/4} SWAP. The three terms commute, so three successive rotations do it.
        swap(ψ) = run_circuit(ψ, ["XX", "YY", "ZZ"], [π/4, π/4, π/4])
        g = cis(π/4)
        for (ones, swapped) in ((Int[], Int[]), ([1], [2]), ([2], [1]), ([1, 2], [1, 2]))
            ψ = swap(basis_state(2, ones))
            @test close(ψ, g * basis_state(2, swapped))
            z1 = 1 in swapped ? -1 : 1
            z2 = 2 in swapped ? -1 : 1
            @test isapprox(expect(ψ, "ZI"), z1; atol = ATOL)
            @test isapprox(expect(ψ, "IZ"), z2; atol = ATOL)
        end
        # |a⟩⊗|b⟩ with |a⟩ = exp(i·0.3·X)|0⟩ on qubit 1, Bloch (0, sin 0.6, cos 0.6),
        # and |b⟩ = exp(i·0.7·Y)|0⟩ on qubit 2, Bloch (−sin 1.4, 0, cos 1.4).
        ba = [0, sin(0.6), cos(0.6)]
        bb = [-sin(1.4), 0, cos(1.4)]
        ψ = run_circuit(zero_state(2), ["XI", "IY"], [0.3, 0.7])
        @test close(collect(bloch(ψ, 1)), ba)
        @test close(collect(bloch(ψ, 2)), bb)
        φ = swap(ψ)
        @test close(collect(bloch(φ, 1)), bb)
        @test close(collect(bloch(φ, 2)), ba)
        # The swapped state is the product |b⟩⊗|a⟩ built directly, times e^{iπ/4}.
        @test close(φ, g * run_circuit(zero_state(2), ["YI", "IX"], [0.7, 0.3]))
        # Product state: ⟨P⊗Q⟩ = ⟨P⟩⟨Q⟩, e.g. ⟨Y₁X₂⟩ before = sin 0.6 · (−sin 1.4),
        # and after the swap ⟨X₁Y₂⟩ takes that value.
        @test isapprox(expect(ψ, "YX"), -sin(0.6) * sin(1.4); atol = ATOL)
        @test isapprox(expect(φ, "XY"), -sin(0.6) * sin(1.4); atol = ATOL)
    end

    @testset "order of gates" begin
        # gates[1] = (X, a) acts first, then (Z, b): ψ = U_Z U_X |0⟩ with
        # U_X = [[cos a, i sin a], [i sin a, cos a]], U_Z = diag(e^{ib}, e^{−ib}),
        # so ψ = (e^{ib} cos a, i e^{−ib} sin a).
        a, b = 0.4, 0.9
        UX = ComplexF64[cos(a) im*sin(a); im*sin(a) cos(a)]
        UZ = ComplexF64[cis(b) 0; 0 cis(-b)]
        e0 = ComplexF64[1, 0]
        ψ = run_circuit(zero_state(1), ["X", "Z"], [a, b])
        @test close(ψ, matvec(UZ, matvec(UX, e0)))
        @test close(ψ, [cis(b) * cos(a), im * cis(-b) * sin(a)])
        # The other order gives U_X U_Z|0⟩ = e^{ib}(cos a, i sin a), a different state.
        χ = run_circuit(zero_state(1), ["Z", "X"], [b, a])
        @test close(χ, matvec(UX, matvec(UZ, e0)))
        @test close(χ, cis(b) * [cos(a), im * sin(a)])
        @test !close(ψ, χ)
        # Physical difference: ⟨X⟩ = 2Re(conj(ψ₀)ψ₁) = 2 cos a sin a Re(i e^{−2ib})
        # = sin 2a sin 2b for ψ, but 0 for χ.
        @test isapprox(expect(ψ, "X"), sin(2a) * sin(2b); atol = ATOL)
        @test isapprox(expect(χ, "X"), 0; atol = ATOL)
        # run_circuit does not mutate its input and equals successive rotate calls.
        ψ0 = zero_state(1)
        @test run_circuit(ψ0, ["X", "Z"], [a, b]) == rotate(rotate(ψ0, "X", a), "Z", b)
        @test ψ0 == zero_state(1)
        @test run_circuit(ψ0, String[], Float64[]) == ψ0
    end

    @testset "qubit ordering" begin
        ψ = basis_state(3, [2])     # |0 1 0⟩: only qubit 2 is in |1⟩
        @test expect(ψ, "ZII") == 1
        @test expect(ψ, "IZI") == -1
        @test expect(ψ, "IIZ") == 1
        @test bloch(ψ, 1) == (0, 0, 1)
        @test bloch(ψ, 2) == (0, 0, -1)
        @test bloch(ψ, 3) == (0, 0, 1)
        # "XII" flips qubit 1 only: |010⟩ → |110⟩.
        @test apply_pauli(ψ, "XII") == basis_state(3, [1, 2])
        @test apply_pauli(ψ, "IIX") == basis_state(3, [2, 3])
        # Documented layout: qubit 1 is the least significant bit of the index, so the
        # matrix of A₁⊗A₂⊗A₃ is kron(A₃, A₂, A₁) (kron's first factor is most significant).
        v = [complex(sin(k), cos(2k + 1)) for k in 1:8]
        v ./= sqrt(norm2(v))
        mats = Dict('I' => I2, 'X' => X2, 'Y' => Y2, 'Z' => Z2)
        for P in ("XIZ", "ZIX", "YXI", "IYZ", "ZZY", "YYY", "XYZ")
            M = kron(mats[P[3]], mats[P[2]], mats[P[1]])
            Mv = matvec(M, v)
            @test close(apply_pauli(v, P), Mv)
            @test isapprox(expect(v, P), real(sum(conj.(v) .* Mv)); atol = ATOL)
        end
        # "XIZ" is X on qubit 1 and Z on qubit 3, and differs from "ZIX" on this state.
        @test !close(apply_pauli(v, "XIZ"), apply_pauli(v, "ZIX"))
        # bloch(ψ, j) addresses qubit j: rotate only qubit 3 by exp(i·0.2·X).
        φ = rotate(zero_state(3), "IIX", 0.2)
        @test close(collect(bloch(φ, 1)), [0, 0, 1])
        @test close(collect(bloch(φ, 2)), [0, 0, 1])
        @test close(collect(bloch(φ, 3)), [0, sin(0.4), cos(0.4)])
    end

    @testset "Pauli algebra" begin
        v = [complex(cos(3k), sin(k * k)) for k in 1:4]
        ps = all_paulis(2)
        @test length(ps) == 16 && length(unique(ps)) == 16
        @test ps[1] == "II" && ps[2] == "XI" && ps[5] == "IX" && ps[16] == "ZZ"
        @test all_paulis(1) == ["I", "X", "Y", "Z"]
        for P in ps      # P² = 1
            @test close(apply_pauli(apply_pauli(v, P), P), v)
        end
        # XY = iZ, YZ = iX, ZX = iY (the rightmost operator acts first).
        w = ComplexF64[0.6 + 0.1im, -0.3 + 0.7im]
        @test close(apply_pauli(apply_pauli(w, "Y"), "X"), im * apply_pauli(w, "Z"))
        @test close(apply_pauli(apply_pauli(w, "Z"), "Y"), im * apply_pauli(w, "X"))
        @test close(apply_pauli(apply_pauli(w, "X"), "Z"), im * apply_pauli(w, "Y"))
        # Y = iXZ (CONVENTIONS §1).
        @test close(apply_pauli(w, "Y"), im * apply_pauli(apply_pauli(w, "Z"), "X"))
        # Anticommutation: XZ = −ZX.
        @test close(apply_pauli(apply_pauli(w, "Z"), "X"), -apply_pauli(apply_pauli(w, "X"), "Z"))
        # |ψ⟩⟨ψ| = 2^{-N} Σ_P ⟨P⟩ P and Tr(PQ) = 2^N δ_PQ, so Tr ρ² = 1 gives Σ_P ⟨P⟩² = 2^N.
        ψ = run_circuit(zero_state(3), ["XXI", "IYZ", "ZIX", "YIY", "IZI", "XYZ"],
                        [0.3, -0.7, 1.1, 0.45, -0.2, 2.0])
        ex = all_expectations(ψ)
        @test length(ex) == 64
        @test isapprox(sum(abs2, ex), 8; atol = ATOL)
        @test ex == [expect(ψ, P) for P in all_paulis(3)]
        @test isapprox(ex[1], 1; atol = ATOL)       # "III"
        # For |000⟩: ⟨P⟩ = 1 if P contains only I and Z (8 strings), else 0.
        ex0 = all_expectations(zero_state(3))
        @test ex0 == [all(c -> c in "IZ", P) ? 1.0 : 0.0 for P in all_paulis(3)]
    end

    @testset "norm and identity" begin
        ψ = zero_state(4)
        @test expect(ψ, "IIII") == 1
        gates = ["XYZI", "IZZY", "-YIIX", "ZXXZ", "IIYI"]
        angles = [0.37, -1.2, π/4, 2.9, -π/2]
        for (P, ϕ) in zip(gates, angles)
            ψ = rotate(ψ, P, ϕ)
            @test isapprox(norm2(ψ), 1; atol = ATOL)
            @test isapprox(expect(ψ, "IIII"), 1; atol = ATOL)
        end
        @test close(ψ, run_circuit(zero_state(4), gates, angles))
        # Undoing the circuit (reverse order, opposite angles) returns |0000⟩.
        @test close(run_circuit(ψ, reverse(gates), -reverse(angles)), zero_state(4))
        # expect is not renormalised: ⟨2ψ|I|2ψ⟩ = 4, ⟨2·0|Z|2·0⟩ = 4.
        @test isapprox(expect(2ψ, "IIII"), 4; atol = ATOL)
        @test expect(2 * zero_state(1), "Z") == 4
    end

    @testset "error handling" begin
        ψ = zero_state(2)
        for bad in ("XA", "xx", "X", "XXX", "", "+", "X-", "+-XX", "X X", "iXX")
            @test_throws ArgumentError apply_pauli(ψ, bad)
            @test_throws ArgumentError expect(ψ, bad)
            @test_throws ArgumentError rotate(ψ, bad, 0.1)
        end
        @test_throws ArgumentError run_circuit(ψ, ["XX", "XQ"], [0.1, 0.2])
        @test_throws ArgumentError run_circuit(ψ, ["XX"], [0.1, 0.2])
        @test_throws ArgumentError basis_state(2, [3])
        @test_throws ArgumentError basis_state(2, [0])
        @test_throws ArgumentError bloch(ψ, 3)
        @test_throws ArgumentError expect(zeros(ComplexF64, 3), "X")
    end
end
