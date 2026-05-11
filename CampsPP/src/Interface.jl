import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

function stringtopauli(str::String)
    char_to_sym = Dict('I' => :I, 'X' => :X, 'Y' => :Y, 'Z' => :Z)
    nq = length(str)
    paulis = Symbol[]
    qinds = Int[]
    for (i, c) in enumerate(str)
        s = char_to_sym[c]
        if s !== :I
            push!(paulis, s)
            push!(qinds, i)
        end
    end
    return pp.PauliString(nq, paulis, qinds)
end

""" ```rotation_circuit(t, N)```

Generate t N-qubit random Pauli rotation gates e^(iϕP), with phases ϕ returned separately.
"""
function rotation_circuit(t::Integer, Nqubits::Integer)
  randstrings = [pp.inttosymbol(rand(0:4^Nqubits-1), Nqubits) for _ in 1:t]
  return rotation_circuit(randstrings, Nqubits)
end

""" ```rotation_circuit(Ps::Vector{Vector{Symbol}}, N)```

Generate N-qubit Pauli rotation gates e^(iϕP) with the given Ps, \
with random phases ϕ returned separately.
"""
function rotation_circuit(Ps::Vector{<:Vector}, Nqubits::Integer)
  randphases = 2π * rand(Float64, (length(Ps),)) # Exponential phases, ie. -1/2 * rotation angles
  return rotation_circuit(Ps, randphases, Nqubits)
end

""" ```rotation_circuit(ϕs, N)```

Generate N-qubit random Pauli rotation gates e^(iϕP), with given phases ϕ returned separately.
"""
function rotation_circuit(ϕs::Vector{<:Real}, Nqubits::Integer)
  randstrings = [pp.inttosymbol(rand(0:4^Nqubits-1), Nqubits) for _ in 1:length(ϕs)]
  return rotation_circuit(randstrings, ϕs, Nqubits)
end

""" ```rotation_circuit(Ps::Vector{Vector{Symbol}}, ϕs, N)```

Generate N-qubit Pauli rotation gates e^(iϕP), with given phases ϕ returned separately.
"""
function rotation_circuit(Ps::Vector{<:Vector}, ϕs::Vector{<:Real}, Nqubits::Integer)
  qinds = [collect(1:Nqubits) for _ in Ps]
  rotations = pp.PauliRotation.(Ps, qinds)
  return rotations, ϕs
end

""" ```random_paulistring(N)```

Generate a random N-qubit Pauli string.
"""
random_paulistring(N::Integer) = stringtopauli(join(rand(["I","X","Y","Z"], N)))

"""```random_rotation(N)```
Generate a random N-qubit Pauli rotation gate e^(iϕP), with ϕ returned separately."""
function random_rotation(Nqubits)
  gate = PauliOperator(random_paulistring(Nqubits))
  phase = 2π * rand(Float64)
  return gate, phase
end
