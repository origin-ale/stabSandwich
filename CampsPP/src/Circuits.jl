import PauliPropagation as pp
using QuantumClifford
using .DisentangleCAMPS
using Random: AbstractRNG, default_rng, randperm

# == Doping ===================================================================

""" ```doping_mask_count([rng,] phases, n)```

Count-mode doping mask over the XX-YY-ZZ gate triples of an XXZ circuit: a \
```BitVector``` with one entry per triple (T = ```length(phases) ÷ 3```) and exactly n \
```true``` entries, at n distinct triples drawn uniformly without replacement over the \
whole circuit. n must be an integer with 0 ≤ n ≤ T, otherwise ```ArgumentError```. \
The mask depends only on the RNG state (```rng``` defaults to the global RNG). Pass it \
to ```xy_magic```/```z_magic```; n counts triples, so ```:xy``` dopes 2n gates."""
function doping_mask_count(rng::AbstractRNG, phases, n)
  T = length(phases) ÷ 3
  (n isa Integer && !(n isa Bool)) ||
    throw(ArgumentError("doping count must be an integer, got $(repr(n))"))
  0 <= n <= T ||
    throw(ArgumentError("doping count must lie in 0:$T (number of triples), got $n"))
  mask = falses(T)
  mask[randperm(rng, T)[1:n]] .= true
  return mask
end
doping_mask_count(phases, n) = doping_mask_count(default_rng(), phases, n)

""" ```doping_mask([rng,] phases, p)```

Probability-mode doping mask (historical: kept to reproduce p-based runs, e.g. the \
T11 replay of Phase 0; the scripts use ```doping_mask_count```). A ```BitVector``` \
with one entry per XX-YY-ZZ triple, entry m ```true``` when the m-th of \
```rand(rng, T)``` is < p. Pass it to ```xy_magic```/```z_magic```, so that the same \
mask reproduces the same doping regardless of the RNG state at call time."""
function doping_mask(rng::AbstractRNG, phases, p)
  return rand(rng, length(phases) ÷ 3) .< p
end
doping_mask(phases, p) = doping_mask(default_rng(), phases, p)

""" ```xy_magic(phases, mask; magicphase = -π/8)```

Dope an XX-YY-ZZ circuit (see ```xxz_circuit```) with magic: in each triple selected \
by the boolean ```mask``` (see ```doping_mask_count```, ```doping_mask```), replace the \
XX and YY phases by ```magicphase```. Returns a new phase vector; ```phases``` is not \
modified."""
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

""" ```z_magic(phases, mask; magicphase = -π/8)```

As ```xy_magic```, but replace the ZZ phase of each doped triple by ```magicphase```."""
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

Generate t N-qubit random Pauli rotation gates e^(iϕP): uniformly random strings on all N \
qubits (identity included) and phases ϕ uniform in [0, 2π). Returns \
```(rotations, phases)```; ```rng``` defaults to the global RNG."""
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

Generate t layers (cycles) of an N-qubit Floquet XXZ circuit: in each layer, the gate \
exp[i(ϕ(XX+YY)+θZZ)] on every bond in ```bricklayer_qinds(N)``` order, each as the three \
rotations XX, YY, ZZ (they commute) with phases ϕ, ϕ, θ. Returns \
```(rotations, phases)```: a vector of ```pp.PauliRotation```s and the matching phases, \
laid out XX-YY-ZZ XX-YY-ZZ … (the triples that ```xy_magic```/```z_magic``` dope)."""
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

""" ```bricklayer_qinds(N)```

The nearest-neighbour qubit pairs of one brickwork layer on N qubits: first \
[1, 2], [3, 4], …, then [2, 3], [4, 5], …"""
function bricklayer_qinds(N::Integer)
  first_half = [[i, i+1] for i in 1:2:N-1]
  second_half = [[i, i+1] for i in 2:2:N-1]
  bonds = []
  append!(bonds, first_half)
  append!(bonds, second_half)
  return bonds
end
