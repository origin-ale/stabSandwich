using Test
import CliffordMPS as cmps
import PauliPropagation as pp
using PauliPropagation: getpauli
using QuantumClifford

include("../CiffordMPSPauliPropagation.jl")

# Independent encoding reference:
#   PauliPropagation packs qubit i at bits 2*(i-1) and 2*(i-1)+1 with I=0, X=1, Y=2, Z=3.
#   Multi-qubit integer = sum over i of pauli_i * 4^(i-1).
#   QuantumClifford P"ABC" assigns A→qubit1, B→qubit2, C→qubit3.

@testset "CiffordMPSPauliPropagation" begin

    @testset "_inttostabilizer single Pauli string" begin
        # Single-qubit cases: code directly encodes the Pauli
        @test _inttostabilizer(0, 1)[1] == P"I"
        @test _inttostabilizer(1, 1)[1] == P"X"
        @test _inttostabilizer(2, 1)[1] == P"Y"
        @test _inttostabilizer(3, 1)[1] == P"Z"

        # 3-qubit XYZ: X@q1=1, Y@q2=2*4=8, Z@q3=3*16=48 → 57
        stab_xyz = _inttostabilizer(57, 3)
        @test length(stab_xyz) == 1
        @test stab_xyz[1] == P"XYZ"

        # 3-qubit IZX: I@q1=0, Z@q2=3*4=12, X@q3=1*16=16 → 28
        stab_izx = _inttostabilizer(28, 3)
        @test length(stab_izx) == 1
        @test stab_izx[1] == P"IZX"

        # 2-qubit YZ: Y@q1=2, Z@q2=3*4=12 → 14
        stab_yz = _inttostabilizer(14, 2)
        @test length(stab_yz) == 1
        @test stab_yz[1] == P"YZ"
    end

    @testset "_inttostabilizer vector of Pauli strings" begin
        # Two 1-qubit strings: X (code 1) and Z (code 3)
        stab = _inttostabilizer(UInt8[1, 3], 1)
        @test length(stab) == 2
        @test stab[1] == P"X"
        @test stab[2] == P"Z"

        # Two 2-qubit strings: XY (1 + 2*4 = 9) and ZI (3 + 0*4 = 3)
        stab2 = _inttostabilizer(UInt8[9, 3], 2)
        @test length(stab2) == 2
        @test stab2[1] == P"XY"
        @test stab2[2] == P"ZI"

        # Single-element vector
        stab3 = _inttostabilizer(UInt8[2], 1)
        @test length(stab3) == 1
        @test stab3[1] == P"Y"
    end

    @testset "cmps.PauliSum from pp.PauliString" begin
        # 1 qubit, 2.5 * X
        pstr = pp.PauliString(1, [:X], [1], 2.5)
        result = cmps.PauliSum(pstr)
        @test length(result.coeffs) == 1
        @test result.coeffs[1] == 2.5
        @test length(result.ops) == 1
        @test result.ops[1] == P"X"

        # 2 qubits, 1.0 * XZ (default coefficient)
        pstr2 = pp.PauliString(2, [:X, :Z], [1, 2])
        result2 = cmps.PauliSum(pstr2)
        @test length(result2.coeffs) == 1
        @test result2.coeffs[1] == 1.0
        @test result2.ops[1] == P"XZ"

        # 3 qubits, -3.0 * YXZ
        pstr3 = pp.PauliString(3, [:Y, :X, :Z], [1, 2, 3], -3.0)
        result3 = cmps.PauliSum(pstr3)
        @test length(result3.coeffs) == 1
        @test result3.coeffs[1] == -3.0
        @test result3.ops[1] == P"YXZ"

        # Identity string on 2 qubits, coefficient 0.5
        pstr_id = pp.PauliString(2, [:I, :I], [1, 2], 0.5)
        result_id = cmps.PauliSum(pstr_id)
        @test length(result_id.coeffs) == 1
        @test result_id.coeffs[1] == 0.5
        @test result_id.ops[1] == P"II"
    end

    @testset "cmps.PauliSum from pp.PauliSum" begin
        # Single term: 1 qubit, 2.5 * X
        psum_single = pp.PauliSum(1, Dict(UInt8(1) => 2.5))
        result = cmps.PauliSum(psum_single)
        @test length(result.coeffs) == 1
        @test result.coeffs[1] == 2.5
        @test result.ops[1] == P"X"

        # Two terms: 1 qubit, 2.0*X + (-1.0)*Z
        # Dict ordering is not guaranteed, so verify with a lookup
        psum2 = pp.PauliSum(1, Dict(UInt8(1) => 2.0, UInt8(3) => -1.0))
        result2 = cmps.PauliSum(psum2)
        @test length(result2.coeffs) == 2
        coeff_by_op = Dict(result2.ops[i] => result2.coeffs[i] for i in 1:2)
        @test coeff_by_op[P"X"] == 2.0
        @test coeff_by_op[P"Z"] == -1.0

        # Two terms: 2 qubits, 1.0*XY + 0.5*ZI
        # XY: 1 + 2*4 = 9, ZI: 3 + 0*4 = 3
        psum3 = pp.PauliSum(2, Dict(UInt8(9) => 1.0, UInt8(3) => 0.5))
        result3 = cmps.PauliSum(psum3)
        @test length(result3.coeffs) == 2
        coeff_by_op3 = Dict(result3.ops[i] => result3.coeffs[i] for i in 1:2)
        @test coeff_by_op3[P"XY"] == 1.0
        @test coeff_by_op3[P"ZI"] == 0.5
    end

end
