using Revise

import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

# == Doping ===================================================================

""" ```dopeMagic(N, gates, phases, layer_ends, dope_syms, dope_inds, p; [magicphase])```

Dope the N-qubit circuit with magic rotations by adding a magic rotation, \
default π/8, with the given symbols and indices after each layer with \
probability p. `dope_inds` may be a callable, in which case it is invoked \
once per inserted rotation to determine the indices (allowing per-layer \
randomization)."""
function dopeMagic(N, gates, phases, layer_ends, dope_syms, dope_inds, p; magicphase = π/8)
  newgates = copy(gates)
  newphases = copy(phases)
  newends = copy(layer_ends)
  for i in eachindex(layer_ends)
    if rand() < p
      inds = dope_inds isa Function ? dope_inds() : dope_inds
      insert!(newgates, newends[i]+1, pp.PauliRotation(dope_syms, inds))
      insert!(newphases, newends[i]+1, magicphase)
      newends[i:end] .+= 1
    end
  end
  return newgates, newphases, newends
end

""" ```subMagic(phases, p; [magicphase])```

Dope a circuit with magic by substituting some phases with magic ones,\
default π/8, with probability p."""
function subMagic(phases, p; magicphase = π/8)
  newphases = copy(phases)
  for i in eachindex(phases)
    if rand() < p
      newphases[i] = magicphase
    end
  end
  return newphases
end

""" ```dopeT(N, gates, phases, layer_ends, p)```

Dope the N-qubit circuit with T gates by adding one on a random index \
after each layer with probability p."""
dopeT(N, gates, phases, layer_ends, p) =
  dopeMagic(N, gates, phases, layer_ends, [:Z], () -> rand(1:N), p)

# == Random rotations ===================================================================

""" ```rotation_circuit(t, N)```

Generate t N-qubit random Pauli rotation gates e^(iϕP), with phases ϕ returned separately."""
function rotation_circuit(t::Integer, N::Integer)
  randstrings = [random_paulistr_sym(N) for _ in 1:t]
  return rotation_circuit(randstrings, N)
end

""" ```rotation_circuit(Ps::Vector{Vector{Symbol}}, N)```

Generate N-qubit Pauli rotation gates e^(iϕP) with the given Ps, \
with random phases ϕ returned separately."""
function rotation_circuit(Ps::Vector{<:Vector}, N::Integer)
  t = length(Ps)
  randphases = 2π * rand(Float64, (t,)) # Exponential phases, ie. -1/2 * rotation angles
  return rotation_circuit(Ps, randphases, N)
end

""" ```rotation_circuit(ϕs, N)```

Generate N-qubit random Pauli rotation gates e^(iϕP), with given phases ϕ returned separately."""
function rotation_circuit(ϕs::Vector{<:Real}, N::Integer)
  randstrings = [random_paulistr_sym(N) for _ in 1:t]
  return rotation_circuit(randstrings, ϕs, N)
end

""" ```rotation_circuit(Ps::Vector{Vector{Symbol}}, ϕs, N)```

Generate N-qubit Pauli rotation gates e^(iϕP), with given phases ϕ returned separately."""
function rotation_circuit(Ps::Vector{<:Vector}, ϕs::Vector{<:Real}, N::Integer)
  qinds = [collect(1:N) for _ in Ps]
  rotations = pp.PauliRotation.(Ps, qinds)
  return rotations, ϕs
end

""" ```random_paulistr_sym(N)```

Generate a random N-qubit Pauli string as a vector of N symbols \
```:I```, ```:X```, ```:Y```, ```:Z``` \
(compatible with PauliPropagation.jl)."""
random_paulistr_sym(N::Integer) = pp.inttosymbol(rand(0:BigInt(4)^N-1), N) 

# == XXZ ================================================================================

""" ```xxz_circuit(ϕ, θ, t, N)```

Generate t layers of a N-qubit Floquet-trotterized XXZ \
(ie. exp[i(ϕ(XX+YY)+θ(ZZ))]) circuit, laid out as\
XX-YY-ZZ XX-YY-ZZ …"""
function xxz_circuit(ϕ, θ, t, N)
  rots = pp.PauliRotation[]
  phases = Real[]
  qinds = staircase_qinds(N)
  for c in 1:t
    for i in qinds
    rots_i, phases_i = xxz(ϕ, θ, i)
    append!(rots, rots_i)
    append!(phases, phases_i)
    end
  end
  return rots, phases
end

""" ```xxz(ϕ, θ, qinds)```

Return rotations and phases corresponding to the gate \
exp[i(ϕ(XX+YY)+θ(ZZ))] acting on qinds."""
function xxz(ϕ, θ, qinds)
  rots = [pp.PauliRotation(p, qinds) for p in 
    [[:X, :X], [:Y, :Y], [:Z, :Z]]
  ]
  phases = [ϕ, ϕ, θ]
  return rots, phases
end

""" ```xxz_layer_circuit(ϕ, θ, t, N)```

Generate t layers of a N-qubit Floquet-trotterized XXZ \
(ie. exp[i(ϕ(XX+YY)+θ(ZZ))]) circuit, laid out as\
XX-XX-… YY-YY-… ZZ-ZZ-… …"""
function xxz_layer_circuit(ϕ, θ, t, N)
  rots = pp.PauliRotation[]
  phases = Real[]
  qinds = bricklayer_qinds(N)
  for c in 1:t
    for g in [([:X,:X], ϕ), ([:Y,:Y], ϕ), ([:Z,:Z], θ)]
      for i in qinds
        rot_i = pp.PauliRotation(g[1], i)
        phase_i = g[2]
        push!(rots, rot_i)
        push!(phases, phase_i)
      end
    end
  end
  return rots, phases
end

# == fSim ===============================================================================

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

# == Utility ===================================================================

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

""" ```staircase_qinds(N)````

Return a vector of pairs of indices corresponding to a staircase \
(1-2,2-3,…) N-qubit circuit structure."""
function staircase_qinds(N::Integer)
  return [[i, i+1] for i in 1:N-1]
end

# """ ```triangle_qinds(N)````

# Return a vector of pairs of indices corresponding to a triangle \
# (1-2,3-4,2-3,5-6,7-8,6-7,…) N-qubit circuit structure."""
# function triangle_qinds(N::Integer)
#   bonds = []
#   for i = 1:4:N-3
#     new = [[i,i+1],[i+2,i+3],[i+1,i+2]]
#     append!(bonds,new)
#   end
#   return bonds
# end