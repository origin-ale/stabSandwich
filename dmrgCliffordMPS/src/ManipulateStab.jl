using QuantumClifford

import Base
import Kronecker

# -- Clifford Operator products -------------------------------------------------------------------

function Base.:^(p::PauliOperator, op::CliffordOperator)
    return apply!(Stabilizer([copy(p)]), op)[1]
end

function QuantumClifford.apply!(C::CliffordOperator,
    op::CliffordOperator, sites::AbstractArray{Int,1})
    apply!(Stabilizer(C.tab), op, sites)
    return nothing
end

function Base.:*(oper::Tuple{CliffordOperator,AbstractArray{Int,1}}, C::CliffordOperator)
    op = oper[1]
    sites = oper[2]
    CCp = copy(C)
    QuantumClifford.apply!(CCp, op, sites)
    return CCp
end

function Base.:*(C::CliffordOperator, oper::Tuple{CliffordOperator,AbstractArray{Int,1}})
    op = oper[1]
    sites = oper[2]
    res = one(CliffordOperator, length(C[1]))
    QuantumClifford.apply!(Stabilizer(res.tab), op, sites)
    QuantumClifford.apply!(Stabilizer(res.tab), C)
    return res
end


# -- PauliOperator --------------------------------------------------------------------------------

function QuantumClifford.PauliOperator(p::PauliOperator, N::Int, sites::AbstractArray{Int,1})
    x = fill(false,N)
    z = fill(false,N)
    phase = p.phase[]
    for i in 1:length(sites)
        x[sites[i]] = p[i][1]
        z[sites[i]] = p[i][2]
    end
    return PauliOperator(phase, x, z)
end

# -- 2Qubit density matrix ------------------------------------------------------------------------

sigmas = [
    [1   0  ; 0   1  ],
    [0   1  ; 1   0  ],
    [0  -im ; im  0  ],
    [1   0  ; 0  -1  ]
]

function density_matrix(psi::Stabilizer, sites::AbstractArray{Int,1})
    L = length(sites)
    N = length(psi[1])
    m = fill(0.0, 2^L, 2^L)
    for ind in CartesianIndices(tuple(repeat([1:4],L)...))
        p = fill(0,N)
        p[sites] = [ind[i]-1 for i in 1:L]
        _, _, val = project!(copy(psi), PauliOperator(p))
        if val !== nothing
            m += im^val * Kronecker.kronecker([sigmas[ind[i]] for i in 1:L]...)
        end
    end
    return m/2.0^L
end

# -- apply gates ----------------------------------------------------------------------------------

function applyGate!(psi::Stabilizer, op::CliffordOperator, sites::AbstractArray{Int,1})
    QuantumClifford.apply!(psi, op, sites)
    return nothing
end
applyGate!(psi::Stabilizer, op::CliffordOperator, n::Int) = QuantumClifford.apply!(psi, op, [n])
applyGate!(psi::Stabilizer, op::CliffordOperator, nL::Int, nR::Int) = QuantumClifford.apply!(psi, op, [nL, nR])


function applyGate!(C::CliffordOperator, op::CliffordOperator, sites::AbstractArray{Int,1})
    QuantumClifford.apply!(C, op, sites)
    return nothing
end
applyGate!(C::CliffordOperator, op::CliffordOperator, n::Int) = QuantumClifford.apply!(C, op, [n])
applyGate!(C::CliffordOperator, op::CliffordOperator, n1::Int, n2::Int) = QuantumClifford.apply!(C, op, [n1, n2])


