import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

using Revise

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
  randstrings = [pp.inttosymbol(rand(0:BigInt(4)^Nqubits-1), Nqubits) for _ in 1:t]
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
function random_rotation(Nqubits, ::PauliOperator)
  gate = PauliOperator(random_paulistring(Nqubits))
  phase = 2π * rand(Float64)
  return gate, phase
end

function random_rotation(Nqubits, ::pp.PauliRotation)
  qinds = collect(1:Nqubits)
  pstr = pp.inttosymbol(rand(0:4^Nqubits-1), Nqubits)
  gate = pp.PauliRotation(pstr, qinds)
  angle = 2π * rand(Float64)
  return gate, angle
end

""" ```xxz_circuit(t, N)```

Generate t layers of a N-qubit Floquet-trotterized XXZ circuit with ϕ=π/4 for each gate."""
function xxz_circuit(t, Nqubits)
  rots = pp.heisenbergtrottercircuit(Nqubits, t)
  ϕs = [π/4 for r in rots]
  return rots, ϕs
end

""" ```dopeT!(gates, phases, p)````
Dope the N-qubit circuit with T gates by adding one on a random index \
after each gate with probability p.\
Also return the positions of T gates."""
function dopeT(N, gates, phases, p)
  newgates = copy(gates)
  newphases = copy(phases)
  magic_pos = []
  os = 0
  for i in eachindex(gates)
    if rand() < p
      insert!(newgates, i+os, pp.PauliRotation([:Z], rand(1:N)))
      insert!(newphases, i+os, π/8)
      push!(magic_pos, i+os)
      os += 1
    end
  end
  return newgates, newphases,magic_pos
end

""" ```fSim(θ, ϕ, qinds)```

Return rotations and phases corresponding to the gate \
fSim(θ, ϕ) (Rosenberg 2024 supplementary material) acting on qinds."""
function fSim(θ, ϕ, qinds)
  rots = [pp.PauliRotation(p, qinds) for p in 
    [[:X, :X], [:Y, :Y], [:I, :I], [:Z, :Z], [:I, :Z], [:Z, :I]]
  ]
  phases = [θ/2, θ/2, ϕ/4, ϕ/4, -ϕ/4, -ϕ/4]
  return rots, phases
end

""" ```bricklayer_qinds(N)````

Return a vector of pairs of indices corresponding to a bricklayer \
(even-odd bond layers) N-qubit circuit structure."""
function bricklayer_qinds(N::Integer)
  even_bonds = [[i, i+1] for i in 1:2:N-1]
  odd_bonds = [[i, i+1] for i in 2:2:N-1]
  bonds = []
  append!(bonds, even_bonds)
  append!(bonds, odd_bonds)
  return bonds
end

""" ```fSim_circuit(θ, ϕ, t, N)```

Return rotations and phases for a t-cycle fSim Heisenberg-Trotter-like \
N-qubit circuit (Rosenberg 2024) with the given θ and ϕ."""
function fSim_circuit(θ, ϕ, t::Integer, N::Integer)
  rots = []
  phases = []
  qinds = bricklayer_qinds(N)
  for c in 1:t
    for i in qinds
    rots_i, phases_i = fSim(θ, ϕ, i)
    append!(rots, rots_i)
    append!(phases, phases_i)
    end
  end
  return rots, phases
end