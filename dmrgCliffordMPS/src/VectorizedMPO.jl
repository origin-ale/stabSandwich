using ITensors
using ITensorMPS
using QuantumClifford

import Base

# -- VectorizedMPDO -------------------------------------------------------------------------------

"""Lightweight container for a vectorized MPDO stored as an MPS with an even number of sites."""
mutable struct VectorizedMPDO
    mps::Any

    function VectorizedMPDO(mps)
        if isodd(length(mps))
            throw(ArgumentError("`mps` length must be even, got $(length(mps))"))
        end
        new(mps)
    end
end

function Base.copy(vmpdo::VectorizedMPDO)
    return VectorizedMPDO(copy(vmpdo.mps))
end

function Base.length(vmpdo::VectorizedMPDO)
    return div(length(vmpdo.mps),2)
end
function LinearAlgebra.mul!(vmpdo::VectorizedMPDO, factor::Number)
    vmpdo.mps[1] .*= factor
    return nothing
end

function LinearAlgebra.normalize!(vmpdo::VectorizedMPDO)
    LinearAlgebra.normalize!(vmpdo.mps)
    return nothing
end

function LinearAlgebra.norm(vmpdo::VectorizedMPDO)
    return LinearAlgebra.norm(vmpdo.mps)
end

