using ITensors
using ITensorMPS

import Base


# -- show MPS -------------------------------------------------------------------------------------
function Base.show(io::IO, mps::MPS)
    print(io, dim.(linkinds(mps)))
end

# -- apply gates ----------------------------------------------------------------------------------

"""Apply a single-site gate identified by `opID` to site `n` of an MPS."""
function applyGate!(mps::MPS, opID, n::Int)
    sites = siteinds(mps)
    mps[n] = noprime(mps[n] * op(opID,sites[n]))
    return nothing
end

"""Apply a nearest-neighbor gate and split the result back into an MPS."""
function applyGateNN!(mps::MPS, opID, n::Int; leftOC=true, svd_kwargs...)
    sites = siteinds(mps)
    orthogonalize!(mps, n)
    W = noprime(mps[n] * mps[n+1] * op(opID, sites[n],sites[n+1]))
    indsL = uniqueinds(mps[n], mps[n+1])
    U, S, V = svd(W, indsL; svd_kwargs...)
    if leftOC
        mps[n] = U * S
        mps[n+1] = V
    else
        mps[n] = U
        mps[n+1] = S * V
    end
    return diag(S)
end

"""Apply a gate to two possibly separated sites by moving them together with swaps."""
function applyGate!(mps::MPS, opID, nL::Int, nR::Int; leftOC=true, svd_kwargs...)
    if nL+1 == nR
        return applyGateNN!(mps, opID, nL; leftOC=true, svd_kwargs...)
    end
    
    if leftOC
        orthogonalize!(mps, nL)
        applySwaps!(mps, nL, nR-1; svd_kwargs...)
        applyGate!(mps, opID, nR-1, nR; leftOC=true, svd_kwargs...)
        applySwaps!(mps, nR-1, nL; svd_kwargs...)
    else
        orthogonalize!(mps, nR)
        applySwaps!(mps, nR, nL+1; svd_kwargs...)
        applyGate!(mps, opID, nL, nL+1; leftOC=false, svd_kwargs...)
        applySwaps!(mps, nL+1, nR; svd_kwargs...)
    end
    return nothing
end

"""Apply swap gates between two sites and return the singular-value spectra encountered."""
function applySwaps!(mps::MPS, n1::Int, n2::Int; svd_kwargs...)
    orthogonalize!(mps, n1)
    specs = []
    if n1 < n2
        for n ∈ n1:n2-1
            push!(specs, applyGate!(mps, "SWAP", n, n+1; leftOC=false, svd_kwargs...))
        end
    elseif n2 < n1
        for n ∈ reverse(n2:n1-1)
            push!(specs, applyGate!(mps, "SWAP", n, n+1; leftOC=true, svd_kwargs...))
        end
    end
    return specs
end

# -- spectra --------------------------------------------------------------------------------------

"""Return the singular values across bond `b` after orthogonalizing the MPS there."""
function spectrum!(mps::MPS, b::Int; svd_kwargs...)
    orthogonalize!(mps, b)
    U, S, V = svd(
        mps[b], (linkinds(mps, b-1)..., siteinds(mps, b)...); svd_kwargs...)
    return diag(S)
end

"""Return the singular values associated with applying `opID` on bond `n` and `n+1`."""
function spectrum!(mps::MPS, opID, n::Int; svd_kwargs...)
    sites = siteinds(mps)
    orthogonalize!(mps, n)
    W = noprime(mps[n] * mps[n+1] * op(opID, sites[n], sites[n+1]))
    indsL = uniqueinds(mps[n], mps[n+1])
    U, S, V = svd(W, indsL; svd_kwargs...)
    return diag(S)
end

"""Return all bond singular-value spectra of an MPS, updating the canonical form in place."""
function spectra!(mps; svd_kwargs...)
    N = length(mps)
    orthogonalize!(mps,N)
    orthogonalize!(mps,1)

    spec = []
    for n in 1:N-1
        W = noprime(mps[n] * mps[n+1])
        indsL = uniqueinds(mps[n], mps[n+1])
        U, S, V = svd(W, indsL; svd_kwargs...)
        push!(spec, diag(S))
        mps[n] = U
        mps[n+1] = S * V        
    end
    orthogonalize!(mps,1)
    return spec
end

# -- entanglement entropy -------------------------------------------------------------------------

"""Compute the Shannon-like entropy of a probability vector in base 2."""
function eEntropy(ps::AbstractArray)
    ee = 0.0
    for p in ps
        if p > 0
            ee -= p * log2(p)
        end
    end
    return ee
end

"""Compute bipartite entanglement entropy from the singular values at bond `b`."""
function eEntropy!(mps::MPS, b::Int; svd_kwargs...)
    return eEntropy(spectrum!(mps, b; svd_kwargs...).^2)
end

"""Compute bipartite entanglement entropy after applying a gate on bond `n` and `n+1`."""
function eEntropy!(mps::MPS, opID, n::Int; svd_kwargs...)
    return eEntropy(spectrum!(mps::MPS, opID, n::Int; svd_kwargs...).^2)
end

"""Compute the entanglement entropy across every bond of the MPS."""
function eEntropys!(mps::MPS; svd_kwargs...)
    return eEntropy.([s.^2 for s in spectra!(mps; svd_kwargs...)])
end

"""Return the total entanglement entropy across all bonds of the MPS."""
function eEntropy!(mps::MPS; svd_kwargs...)
    return sum(eEntropys!(mps))
end

"""Create a normalized random qubit MPS with optional bond-dimension control."""
function random_qubit_mps(N; linkdims=nothing, type="full")
    # ensure a sensible default for `linkdims` before delegating to random_mps
    if linkdims === nothing
        linkdims = 2^(div(N,2))
    end

    full = nothing
    if type=="ITensors"
        mps = random_mps(ComplexF64, siteinds("Qubit",N); linkdims=linkdims)
        orthogonalize!(mps, N)
        orthogonalize!(mps, 1)
        return mps
    elseif type=="full"
        full = true
    elseif type=="small"
        full = false
    else
        throw(ErrorException("type $type is not implemented."))
    end
    sites = siteinds("Qubit",N)
    linkdims_list = nothing
    if full
        linkdims_list = [linkdims for n∈1:N-1]
    else
        linkdims_list = [min(2^min(n, N-n), linkdims) for n∈1:N-1]
    end
    links = [Index(dim,"link") for dim∈linkdims_list]

    Ts = ITensor[]
    push!(Ts, random_itensor(ComplexF64, sites[1],links[1]))
    for n∈2:N-1
        push!(Ts, random_itensor(ComplexF64, links[n-1],sites[n],links[n]))
    end
    push!(Ts, random_itensor(ComplexF64, links[N-1], sites[N]))

    mps = MPS(Ts)
    normalize!(mps)
    orthogonalize!(mps, N)
    orthogonalize!(mps, 1)
    return mps
end

