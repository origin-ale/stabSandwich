using CliffordMPS
using ITensors
using ITensorMPS
using QuantumClifford

import Base

# -- CAMPS ----------------------------------------------------------------------------------------

"""
`CAMPS` stores a Clifford-augmented MPS representation of a qubit state.

The field `mps` contains the non-Clifford tensor-network part, while `Cdag`
stores the Clifford dressing frame.
"""
mutable struct CAMPS
    # The state is |psi> = C|mps>, i.e. |psi> = inv(Cdag)|mps>
    mps::MPS
    Cdag::CliffordOperator
end

"""Construct a dressed state from an existing MPS and an identity Clifford frame."""
function CAMPS(mps::MPS)
    N = length(mps)
    return CAMPS(mps, one(CliffordOperator,N))
end

"""Construct a dressed state from a Clifford frame and a computational-basis MPS."""
function CAMPS(Cdag::CliffordOperator)
    N = length(Cdag[1])
    return CAMPS(MPS(ComplexF64,siteinds("Qubit",N),"0"), Cdag)
end

"""Construct a dressed state on `N` qubits in the computational basis."""
function CAMPS(N::Int)
    return CAMPS(MPS(ComplexF64,siteinds("Qubit",N),"0"), one(CliffordOperator,N))
end

function Base.copy(camps::CAMPS)
    return CAMPS(copy(camps.mps),copy(camps.Cdag))
end

function Base.length(camps::CAMPS)
    return length(camps.mps)
end

function Base.show(io::IO, camps::CAMPS)
    println(io, "CAMPS for ", length(camps), " Qubits.")
    println(io, "----- Cdag -----")
    println(io, camps.Cdag)
    println(io, "----- MPS  -----")
    println(io, camps.mps)
end


function Base.:≈(campsA::CAMPS, campsB::CAMPS)
    @assert length(campsA) == length(campsB)
    for n in 1:length(campsA)
        for p in [P"X", P"Y", P"Z"]
            if 2.0+expectation(campsA, p, [n]) ≉ 2.0+expectation(campsB, p, [n])
                return false
            end
        end
    end
    return true
end

function Base.:≉(campsA::CAMPS, campsB::CAMPS)
    return !(campsA≈campsB)
end

function random_camps(N::Int; debth=nothing, linkdims=nothing, type="ITensors")
    if debth === nothing
        debth = N^2
    end
    if linkdims === nothing
        linkdims = 2^div(N,2)
    end
    
    Cdag = random_clifford_circuit(N, debth)
    mps = random_qubit_mps(N; linkdims=linkdims, type=type)
  
    return CAMPS(mps, Cdag)
end

function LinearAlgebra.mul!(camps::CAMPS, factor::Number)
    camps.mps[1] .*= factor
    return nothing
end

function LinearAlgebra.normalize!(camps::CAMPS)
    LinearAlgebra.normalize!(camps.mps)
    return nothing
end

function LinearAlgebra.norm(camps::CAMPS)
    return LinearAlgebra.norm(camps.mps)
end

# -- Manipulate CAMPS -----------------------------------------------------------------------------

"""Apply an indexed Clifford gate to both the Clifford frame and the MPS part."""
function transform!(camps::CAMPS, ind::Clifford1Index, n::Int)
    applyGate!(camps.Cdag, ind, n)
    applyGate!(camps.mps, ind, n)
    return nothing
end

"""Apply a two-site indexed Clifford gate to both parts of a `CAMPS` state."""
function transform!(camps::CAMPS, ind::Clifford2Index, nL::Int, nR::Int; svd_kwargs...)
    applyGate!(camps.Cdag, ind, nL, nR)
    applyGate!(camps.mps, ind, nL, nR; svd_kwargs...)
    return nothing
end

# -- Use CAMPS ------------------------------------------------------------------------------------

"""Evaluate a Pauli-sum observable on a dressed state."""
function expectation(camps::CAMPS, h::PauliSum)
    hT = h^camps.Cdag
    return ITensorMPS.expect(camps.mps, hT)
end

"""Evaluate a single Pauli operator on a dressed state."""
function expectation(camps::CAMPS, p::PauliOperator)
    @assert length(p) == length(camps)
    h = PauliSum(Stabilizer([p]))
    return expectation(camps, h)
end

"""Embed a local Pauli operator on `sites` and evaluate it on the dressed state."""
function expectation(camps::CAMPS, p::PauliOperator, sites::AbstractArray{Int,1})
    pFull = PauliOperator(p, length(camps), sites)
    return expectation(camps, pFull)
end

"""Apply a Pauli-sum operator to the dressed state and update the MPS part."""
function applyGate!(camps::CAMPS, h::PauliSum; svd_kwargs...)
    hT = h^camps.Cdag
    camps.mps = ITensorMPS.apply(camps.mps, hT; svd_kwargs...)
    return nothing
end 

"""Apply a Clifford gate to the dressed state frame."""
function applyGate!(camps::CAMPS, op::CliffordOperator,  sites::AbstractArray{Int,1})
    camps.Cdag = camps.Cdag * (inv(op), sites)
    return nothing
end 

"""Apply an indexed two-site Clifford gate to the dressed state frame."""
function applyGate!(camps::CAMPS, ind::Clifford2Index,  sites::AbstractArray{Int,1})
    return applyGate!(camps, clifford_gate(ind), sites)
end 


"""Return the reduced density matrix of the dressed state on `sites`."""
function density_matrix(camps::CAMPS, sites::AbstractArray{Int,1})
    L = length(sites)
    N = length(camps)
    m = fill(0.0, 2^L, 2^L)
    for ind in CartesianIndices(tuple(repeat([1:4],L)...))
        p = [ind[i]-1 for i in 1:L]
        exp = expectation(camps, PauliOperator(p), sites)
        if exp != 0.0
            m += exp * Kronecker.kronecker([sigmas[ind[i]] for i in 1:L]...)
        end
    end
    return m/2.0^L
end

"""Return the stabilizer state represented by the Clifford frame in `camps`."""
function QuantumClifford.Stabilizer(camps::CAMPS)
    return QuantumClifford.canonicalize_clip!(inv(camps.Cdag) * one(Stabilizer,length(camps)))
end

function sEntropy(camps::CAMPS, args...; kwargs...)
"""Forward stabilizer-magic entropy evaluation to the MPS part of `camps`."""
    return sEntropy(camps.mps, args...; kwargs...)
end


# -- ITensor CAMPS extensions ---------------------------------------------------------------------

"""Run DMRG on the Hamiltonian dressed by the current Clifford frame."""
function dmrg!(camps::CAMPS, h::PauliSum; dmrg_kwargs...)
    H = MPO(OpSum(h^camps.Cdag), siteinds(camps.mps))
    energy, camps.mps = ITensorMPS.dmrg(H, camps.mps; dmrg_kwargs...)
    return energy
end

"""Run TDVP on the Hamiltonian dressed by the current Clifford frame."""
function tdvp!(camps::CAMPS, h::PauliSum, t::Number; tdvp_kwargs...)
    H = MPO(OpSum(h^camps.Cdag), siteinds(camps.mps))
    camps.mps = tdvp(H, t, camps.mps; tdvp_kwargs...)
end
