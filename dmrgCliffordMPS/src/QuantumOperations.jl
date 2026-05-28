using CliffordMPS
using QuantumClifford
using LinearAlgebra
using Random

abstract type QuantumOperation end
    
"""Typed token describing a high-level quantum operation understood by `applyQOp!`."""
struct QOp{x} <: QuantumOperation
end

QOp(x) = QOp{x}()

# -- Pauli Sum 1gates -----------------------------------------------------------------------------

QOpPauliSum1Gates = Union{QOp{:SqT}, QOp{:T}, QOp{:P0}, QOp{:P1}}

function CliffordMPS.PauliSum(::QOp{:T}, N::Int, site::Int)
    coeffs = [real(sqrt(sqrt(im))),imag(sqrt(sqrt(im)))*(-1.0im) ]
    ops = Stabilizer([PauliOperator(P"I", N, [site]), PauliOperator(P"Z", N, [site])])
    return CliffordMPS.PauliSum(coeffs, ops)
end

function CliffordMPS.PauliSum(::QOp{:SqT}, N::Int, site::Int)
    coeffs = [real(sqrt(sqrt(sqrt(im)))),imag(sqrt(sqrt(sqrt(im))))*(-1.0im) ]
    ops = Stabilizer([PauliOperator(P"I", N, [site]), PauliOperator(P"Z", N, [site])])
    return CliffordMPS.PauliSum(coeffs, ops)
end

function CliffordMPS.PauliSum(::QOp{:P0}, N::Int, site::Int)
    coeffs = [0.5, 0.5]
    ops = Stabilizer([PauliOperator(P"I", N, [site]), PauliOperator(P"Z", N, [site])])
    return CliffordMPS.PauliSum(coeffs, ops)
end

function CliffordMPS.PauliSum(::QOp{:P1}, N::Int, site::Int)
    coeffs = [0.5, -0.5]
    ops = Stabilizer([PauliOperator(P"I", N, [site]), PauliOperator(P"Z", N, [site])])
    return CliffordMPS.PauliSum(coeffs, ops)
end

"""Apply a high-level Pauli-sum operation to a single site of a `CAMPS` state."""
function applyQOp!(camps::CAMPS, qop::QOpPauliSum1Gates, site::Int; svd_kwargs...)
    ps = PauliSum(qop, length(camps), site)
    return applyGate!(camps, ps; svd_kwargs...)
end


# -- Clifford 1gates -----------------------------------------------------------------------------

QOpClifford1Gates = Union{QOp{:H}, QOp{:S}, QOp{:X}, QOp{:Y}, QOp{:Z}}

QuantumClifford.CliffordOperator(::QOp{:H}) = tHadamard
QuantumClifford.CliffordOperator(::QOp{:S}) = tPhase
QuantumClifford.CliffordOperator(::QOp{:X}) = CliffordOperator(sX)
QuantumClifford.CliffordOperator(::QOp{:Y}) = CliffordOperator(sY)
QuantumClifford.CliffordOperator(::QOp{:Z}) = CliffordOperator(sZ)

"""Apply a high-level one-qubit Clifford operation to a `CAMPS` state."""
function applyQOp!(camps::CAMPS, qop::QOpClifford1Gates, site::Int; kwargs...)
    return applyGate!(camps, CliffordOperator(qop), [site])
end


# -- Clifford 2gates -----------------------------------------------------------------------------

QOpClifford2Gates = Union{QOp{:CNOT}, QOp{:SWAP}}

QuantumClifford.CliffordOperator(::QOp{:CNOT}) = tCNOT
QuantumClifford.CliffordOperator(::QOp{:SWAP}) = tSWAP

"""Apply a high-level two-qubit Clifford operation to a `CAMPS` state."""
function applyQOp!(camps::CAMPS, qop::QOpClifford2Gates, site1::Int, site2::Int; kwargs...)
    return applyGate!(camps, CliffordOperator(qop), [site1, site2])
end


# -- Normalization --------------------------------------------------------------------------------

"""Normalize the MPS part of a `CAMPS` state."""
function applyQOp!(camps::CAMPS, ::QOp{:norm}; kwargs...)
    return normalize!(camps)
end

"""Scale the MPS part of a `CAMPS` state by a scalar factor."""
function applyQOp!(camps::CAMPS, ::QOp{:mul}, factor::Number; kwargs...)
    return mul!(camps, factor)
end

# -- Measurement ----------------------------------------------------------------------------------

"""Projectively measure `site` in the Z basis and renormalize the state."""
function applyQOp!(camps::CAMPS, ::QOp{:measurement}, site::Int; kwargs...)
    exp = real(expectation(camps, P"Z", [site]))
    p = 0.5+0.5*exp
    if rand()<p
        applyQOp!(camps, QOp(:P0), site; kwargs...)
        applyQOp!(camps, QOp(:mul), sqrt(1/p))
        ret = 0
    else
        applyQOp!(camps, QOp(:P1), site; kwargs...)
        applyQOp!(camps, QOp(:mul), sqrt(1/(1-p)))
        ret = 1
    end
    return ret
end

# -- Gate Layers ----------------------------------------------------------------------------------

"""Apply a brickwork layer of the same Clifford gate on alternating nearest-neighbor bonds."""
function applyQOp!(camps::CAMPS, ::QOp{:cliffBrickWork}, ind::Clifford2Index, odd::Bool; kwargs...)
    N=length(camps)
    sites = odd ? (1:2:N-1) : (2:2:N-1)
    
    gate = clifford_gate(ind)
    for site in sites
        applyGate!(camps, gate, [site, site+1])
    end
    return nothing
end

"""Apply a random brickwork layer of Clifford gates on alternating nearest-neighbor bonds."""
function applyQOp!(camps::CAMPS, ::QOp{:randCliffBrickWork}, odd::Bool;
    clifford_set=Clifford2IndexSet(:all)::Clifford2IndexSets, kwargs...)
    N=length(camps)
    
    sites = odd ? (1:2:N-1) : (2:2:N-1)
    
    for site in sites
        num = rand(1:length(clifford_set))-1
        ind = Clifford2Index(clifford_set, num)
        gate = clifford_gate(ind)
        applyGate!(camps, gate, [site, site+1])
    end
    return nothing
end

"""Apply multiple random brickwork layers of Clifford gates."""
function applyQOp!(camps::CAMPS, ::QOp{:randCliffBrickWork}; debth=nothing,
    clifford_set=Clifford2IndexSet(:all)::Clifford2IndexSets, kwargs...)
    N=length(camps)
    
    if debth === nothing
        debth = 2*N
    end

    for t in 1:debth
        applyQOp!(camps, QOp(:randCliffBrickWork), t%2 == 1;
            clifford_set=clifford_set, kwargs...)
    end

    return nothing
end

"""Apply a random Clifford circuit to all sites of a `CAMPS` state."""
function applyQOp!(camps::CAMPS, ::QOp{:randCliffCircuit}; debth=nothing,
    clifford_set=Clifford2IndexSet(:all)::Clifford2IndexSets, kwargs...)
    N=length(camps)

    if debth === nothing
        debth = 2*N^2
    end

    circ = random_clifford_circuit(N,debth; clifford_set=clifford_set)

    applyGate!(camps, circ, collect(1:N))

    return nothing
end

"""Apply a random Clifford circuit to a selected subset of sites."""
function applyQOp!(camps::CAMPS, ::QOp{:randCliffCircuit}, sites::AbstractArray{Int,1};
    debth=nothing, clifford_set=Clifford2IndexSet(:all)::Clifford2IndexSets, kwargs...)
    N = length(camps)
    numSites = length(sites)

    if debth === nothing
        debth = 2*numSites^2
    end

    circ = random_clifford_circuit(numSites, debth; clifford_set=clifford_set)

    applyGate!(camps, circ, sites)

    return nothing
end

"""Insert random `T` gates with the given density and disentangle after each insertion."""
function applyQOp!(camps::CAMPS, ::QOp{:randTGates}, density::Float64;
    disentangle_args=tuple(DisentangleStrategy(:none)), disentangle_kwargs=Dict(), kwargs...)
    N=length(camps)
    disentangle_results = []
    sites = Int[]
    for n in 1:N
        if rand()<density
            applyQOp!(camps, QOp(:T), n; kwargs...)
            push!(disentangle_results,
                  disentangle!(camps, disentangle_args...; disentangle_kwargs..., kwargs...))
            push!(sites, n)
        end
    end
    return sites, disentangle_results
end

"""Insert random Z measurements with the given density and disentangle after each measurement."""
function applyQOp!(camps::CAMPS, ::QOp{:randMeasurements}, density::Float64;
    disentangle_args=tuple(DisentangleStrategy(:none)), disentangle_kwargs=Dict(), kwargs...)
    N=length(camps)
    disentangle_results = []
    sites = Int[]
    outcomes = Int[]
    for n in 1:N
        if rand()<density
            outcome = applyQOp!(camps, QOp(:measurement), n; kwargs...)
            push!(disentangle_results,
                  disentangle!(camps, disentangle_args...; disentangle_kwargs..., kwargs...))
            push!(sites, n)
            push!(outcomes, outcome)
        end
    end
    return sites, outcomes, disentangle_results
end
