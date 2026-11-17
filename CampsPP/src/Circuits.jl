using Revise

import PauliPropagation as pp
using QuantumClifford
using .DisentangleCAMPS
using Random: AbstractRNG, default_rng

# == Doping ===================================================================

""" ```doping_mask(rng, phases, p)```

Generate a random doping mask over the XX-YY-ZZ gate triples of an XXZ circuit: \
a `BitVector` with one entry per triple, each `true` with probability p. \
Pass it to `xy_magic`/`z_magic` to select which gates to dope, so that the same \
mask reproduces the same doping regardless of the RNG state at call time."""
function doping_mask(rng::AbstractRNG, phases, p)
  return rand(rng, length(phases) ÷ 3) .< p
end
doping_mask(phases, p) = doping_mask(default_rng(), phases, p)

""" ```xy_magic(phases, p; [magicphase])```
    ```xy_magic(phases, mask; [magicphase])```

Dope an XX-YY-ZZ circuit with magic \
by adding a magic phase to the XX and YY phases,\
default -π/8. Either draw doped triples with probability p, or dope exactly the \
triples selected by a boolean `mask` (see `doping_mask`)."""
function xy_magic(rng::AbstractRNG, phases, p; magicphase = -π/8)
  newphases = copy(phases)
  for i = 1:3:length(phases)
    if rand(rng) < p
      newphases[i] = magicphase
      newphases[i+1] = magicphase
    end
  end
  return newphases
end
xy_magic(phases, p; kwargs...) = xy_magic(default_rng(), phases, p; kwargs...)
function xy_magic(phases, mask::AbstractVector{Bool}; magicphase = -π/8)
  newphases = copy(phases)
  for k in eachindex(mask)
    if mask[k]
      newphases[3k-2] = magicphase
      newphases[3k-1] = magicphase
    end
  end
  return newphases
end

""" ```z_magic(phases, p; [magicphase])```
    ```z_magic(phases, mask; [magicphase])```

Dope an XX-YY-ZZ circuit with magic \
by adding a magic phase to the Z phases,\
default -π/8. Either draw doped triples with probability p, or dope exactly the \
triples selected by a boolean `mask` (see `doping_mask`)."""
function z_magic(rng::AbstractRNG, phases, p; magicphase = -π/8)
  newphases = copy(phases)
  for i = 3:3:length(phases)
    if rand(rng) < p
      newphases[i] = magicphase
    end
  end
  return newphases
end
z_magic(phases, p; kwargs...) = z_magic(default_rng(), phases, p; kwargs...)
function z_magic(phases, mask::AbstractVector{Bool}; magicphase = -π/8)
  newphases = copy(phases)
  for k in eachindex(mask)
    if mask[k]
      newphases[3k] = magicphase
    end
  end
  return newphases
end

# == Random rotations ===================================================================

""" ```rotation_circuit(t, N)```
    ```rotation_circuit(rng, t, N)```

Generate t N-qubit random Pauli rotation gates e^(iϕP), with phases ϕ returned separately."""
function rotation_circuit(rng::AbstractRNG, t::Integer, N::Integer)
  randstrings = [random_paulistr_sym(rng, N) for _ in 1:t]
  return rotation_circuit(rng, randstrings, N)
end
rotation_circuit(t::Integer, N::Integer) = rotation_circuit(default_rng(), t, N)

""" ```rotation_circuit(Ps::Vector{Vector{Symbol}}, N)```
    ```rotation_circuit(rng, Ps::Vector{Vector{Symbol}}, N)```

Generate N-qubit Pauli rotation gates e^(iϕP) with the given Ps, \
with random phases ϕ returned separately."""
function rotation_circuit(rng::AbstractRNG, Ps::Vector{<:Vector}, N::Integer)
  t = length(Ps)
  randphases = 2π * rand(rng, Float64, (t,)) # Exponential phases, ie. -1/2 * rotation angles
  return rotation_circuit(Ps, randphases, N)
end
rotation_circuit(Ps::Vector{<:Vector}, N::Integer) =
  rotation_circuit(default_rng(), Ps, N)

""" ```rotation_circuit(ϕs, N)```
    ```rotation_circuit(rng, ϕs, N)```

Generate N-qubit random Pauli rotation gates e^(iϕP), with given phases ϕ returned separately."""
function rotation_circuit(rng::AbstractRNG, ϕs::Vector{<:Real}, N::Integer)
  randstrings = [random_paulistr_sym(rng, N) for _ in 1:length(ϕs)]
  return rotation_circuit(randstrings, ϕs, N)
end
rotation_circuit(ϕs::Vector{<:Real}, N::Integer) =
  rotation_circuit(default_rng(), ϕs, N)

""" ```rotation_circuit(Ps::Vector{Vector{Symbol}}, ϕs, N)```

Generate N-qubit Pauli rotation gates e^(iϕP), with given phases ϕ returned separately."""
function rotation_circuit(Ps::Vector{<:Vector}, ϕs::Vector{<:Real}, N::Integer)
  qinds = [collect(1:N) for _ in Ps]
  rotations = pp.PauliRotation.(Ps, qinds)
  return rotations, ϕs
end

""" ```random_paulistr_sym(N)```
    ```random_paulistr_sym(rng, N)```

Generate a random N-qubit Pauli string as a vector of N symbols \
```:I```, ```:X```, ```:Y```, ```:Z``` \
(compatible with PauliPropagation.jl)."""
random_paulistr_sym(rng::AbstractRNG, N::Integer) =
  pp.inttosymbol(rand(rng, 0:BigInt(4)^N-1), N)
random_paulistr_sym(N::Integer) = random_paulistr_sym(default_rng(), N)

# == XXZ ================================================================================

""" ```xxz_circuit(ϕ, θ, t, N)```

Generate t layers of a N-qubit Floquet-trotterized XXZ \
(ie. exp[i(ϕ(XX+YY)+θ(ZZ))]) circuit, laid out as\
XX-YY-ZZ XX-YY-ZZ …"""
function xxz_circuit(ϕ, θ, t, N)
  rots = pp.PauliRotation[]
  phases = Real[]
  qinds = bricklayer_qinds(N)
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
