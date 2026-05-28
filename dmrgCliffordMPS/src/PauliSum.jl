using QuantumClifford
using ITensors
using ITensorMPS
using Printf

import Base

# -- Pauli Operators ------------------------------------------------------------------------------

sigmas = [
    ComplexF64[ 1.0 0.0 ; 0.0 1.0],
    ComplexF64[ 0.0 1.0 ; 1.0 0.0],
    ComplexF64[ 0.0 -1.0im ; 1.0im 0.0],
    ComplexF64[ 1.0 0.0 ; 0.0 -1.0]
]

function QuantumClifford.PauliOperator(p::AbstractVector{Int}; phase=0x0)
    x = map(p -> _x_from_int(p),p)
    z = map(p -> _z_from_int(p),p)
    return PauliOperator(phase, x, z)
end

function Base.:*(op::CliffordOperator, p::PauliOperator)
    s = op * Stabilizer([p])
    return s[1]
end

function apply(p::PauliOperator, op::CliffordOperator, indices)
    s = QuantumClifford.apply!(Stabilizer([p]), op, indices)
    return s[1]
end

function _x_from_int(i::Int)
    return i == 1 || i == 2
end

function _z_from_int(i::Int)
    return i == 2 || i == 3
end

function _int_from_xz(x::Bool, z::Bool)
    i = nothing
    if !x && !z
        i = 0
    elseif x && !z
        i = 1
    elseif x && z
        i = 2
    elseif !x && z
        i = 3
    end
    return i
end


# -- Pauli Hamiltonian ----------------------------------------------------------------------------

"""
Sparse Hamiltonian-like sum of Pauli strings with coefficients.

`PauliSum` stores a coefficient vector and a `Stabilizer` basis of Pauli
operators. It is designed for Clifford conjugation, expectation values, and
conversion to `OpSum` or `MPO` representations.
"""
mutable struct PauliSum
    coeffs::AbstractVector{Number}
    ops::Stabilizer
    function PauliSum(coeffs, ops)
        return new(coeffs, ops)
    end
end

"""Construct the zero Pauli sum on `N` qubits."""
function PauliSum(N::Int)
    return PauliSum([0.0], Stabilizer([PauliOperator(P"_",N,[1])]))
end

"""Construct a Pauli sum with unit coefficients from a stabilizer basis."""
function PauliSum(ops::Stabilizer)
    return PauliSum(fill(1.0,length(ops)), ops)
end

function Base.show(io::IO, h::PauliSum)
    for d in 1:debth(h)
        if typeof(h.coeffs[d]) <: Real
            @printf(io,"%+5e ",h.coeffs[d])
        elseif typeof(h.coeffs[d]) <: Complex
            @printf(io,"(%+4e %+4eim) ",real(h.coeffs[d]), imag(h.coeffs[d]))
        else
            print(io,h.coeffs[d])
        end
        println(io,h.ops[d])
    end
end

function Base.copy(h::PauliSum)
    return PauliSum(copy(h.coeffs),copy(h.ops))
end

"""Return the number of Pauli terms stored in the sum."""
function debth(h::PauliSum)
    return length(h.ops)
end

function Base.length(h::PauliSum)
    return length(h.ops[1])
end

"""Append a `(coefficient, PauliOperator)` pair to the sum."""
function Base.push!(h::PauliSum, coeff_pauli::Tuple{Float64, PauliOperator})
    push!(h.coeffs, coeff_pauli[1])
    h.ops = Stabilizer(
        PauliOperator{Array{UInt8, 0}, Vector{UInt64}}[h.ops[:]..., coeff_pauli[2]])
    return nothing
end

# -- Pauli Hamiltonian -> Stab ---------------------------------------------------------------------

"""Conjugate a Pauli sum by a Clifford operator."""
function Base.:^(h::PauliSum, op::CliffordOperator)
    res = copy(h)
    apply!(res, op)
    return res
end

"""Apply Clifford conjugation to the stored Pauli basis in place."""
function QuantumClifford.apply!(
    h::PauliSum, op::CliffordOperator)
    apply!(h.ops, op)
    return nothing
end

"""Apply Clifford conjugation to the stored Pauli basis on the given sites."""
function QuantumClifford.apply!(
    h::PauliSum, op::CliffordOperator, indices::AbstractArray{Int,1})
    apply!(h.ops, op, indices)
    return nothing
end

"""Evaluate a Pauli-sum observable on a stabilizer state."""
function QuantumClifford.expect(state::Stabilizer, h::PauliSum)
    deb = length(h.coeffs)
    nrg = 0.0
    for j in 1:deb
        # _, _, phase = project!(copy(state), h.ops[j]; keep_result=true)
        # if phase != nothing
        #     nrg += 1.0im^phase * h.coeffs[j]
        # end
        nrg += QuantumClifford.expect(h.ops[j], state) * h.coeffs[j]
    end
    return nrg
end
expectation(state::Stabilizer, h::PauliSum) =  QuantumClifford.expect(state, h)

# -- Pauli Hamiltonian -> MPO ---------------------------------------------------------------------

"""Build the local MPO tensor for site `n` from a Pauli sum."""
function tensor(h::PauliSum, n::Int; coeff=false)
    deb = debth(h)
    N = length(h)
    ten = nothing

    @assert(N>1)

    if n==1
        ten = zeros(ComplexF64,2,2,deb)
    elseif n==N
        ten = zeros(ComplexF64,deb,2,2)
    else
        ten = zeros(ComplexF64,deb,2,2,deb)
    end

    for d in 1:deb
        num = _int_from_xz(h.ops[d][n]...)
        if coeff
            c = h.coeffs[d]
        else
            c = 1.0
        end

        if n==1
            ten[:,:,d] = c * sigmas[1+num]
        elseif n==N
            ten[d,:,:] = c * sigmas[1+num]
        else
            ten[d,:,:,d] = c * sigmas[1+num]
        end
    end

    return ten
end

"""Convert a `PauliSum` to an `ITensorMPS.MPO` on the provided sites."""
function ITensorMPS.MPO(sites::Vector{<:Index}, h::PauliSum; ncoeff=0)
    N = length(h)
    deb = debth(h)
    links = [Index(deb, "Link,l=$ii") for ii in 1:(N - 1)]
    
    @assert(N>1)
    @assert(length(sites)==N)

    tens = []
    push!(tens, ITensor(tensor(h,1;coeff=(ncoeff==1)),sites[1],sites[1]',links[1]))
    for n in 2:N-1
        push!(tens,
            ITensor(tensor(h,n;coeff=(ncoeff==n)),links[n-1], sites[n],sites[n]',links[n]))
    end
    push!(tens, ITensor(tensor(h,N;coeff=(ncoeff==N)),links[N-1],sites[N],sites[N]'))

    return MPO(tens, 0, N)
end

"""Update selected MPO sites after a Pauli-sum coefficient or term change."""
function updateMPO!(mpo::MPO, h::PauliSum, n::Vector{Int}; ncoeff=0)
    N = length(h)
    sites =  vcat(siteinds(mpo; plev=0)...)

    for nn in n
        if nn==1
            mpo[1] = ITensor(tensor(h,1;coeff=(ncoeff==1)),sites[1],sites[1]',linkind(mpo,1))
        elseif nn==N
            mpo[N] = ITensor(tensor(h,N;coeff=(ncoeff==N)),linkind(mpo,N-1), sites[N],sites[N]')
        else
            mpo[nn] = ITensor(tensor(h,nn;coeff=(ncoeff==nn)),
            linkind(mpo,nn-1),sites[nn],sites[nn]',linkind(mpo,nn))
        end
    end
    mpo.llim = minimum([minimum(n .- 1), mpo.llim])
    mpo.rlim = maximum([maximum(n .+ 1), mpo.rlim])
    return mpo
end

"""Convert a `PauliSum` to an `ITensorMPS.OpSum`."""
function ITensorMPS.OpSum(h::PauliSum)
    ops = ITensorMPS.OpSum()
    opNames = ["I", "X", "Y", "Z"]
    N = length(h)
    deb = debth(h)
    for d in 1:deb
        l = Any[]
        push!(l,h.coeffs[d]*(1.0im^h.ops[d].phase[]))
        for n in 1:N
            num = _int_from_xz(h.ops[d][n]...)
            push!(l,opNames[num+1])
            push!(l,n)
        end
        add!(ops, l...)
    end
    return ops
end

"""Evaluate a Pauli sum on an MPS by converting it to an MPO first."""
function ITensorMPS.expect(mps::MPS, h::PauliSum)
    sites = siteinds(mps)
    mpo = ITensorMPS.MPO(OpSum(h),sites)
    return ITensorMPS.inner(mps',mpo,mps)
end
expectation(mps::MPS, h::PauliSum) = ITensorMPS.expect(mps, h)

"""Apply a Pauli-sum MPO to an MPS."""
function ITensorMPS.apply(mps::MPS, h::PauliSum; svd_kwargs...)
    sites = siteinds(mps)
    mpo = ITensorMPS.MPO(OpSum(h),sites)
    return ITensorMPS.apply(mpo, mps; svd_kwargs...)
end

