using Revise
using CampsPP
using .DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps
import MagicTensors as mt

# == States =============================================================================

""" ```domainwallstate([T,] rng, N, μ)```

Return an N-qubit product state of type T from the Rosenberg et al. domain wall \
ensemble with contrast parameter μ, and its `onebitinds`. T is `cmps.CAMPS` (the
default while the port is staged) or `MagicTensors.QubitCAMPS`; the random draws
do not depend on T.
"""
domainwallstate(rng, N, μ) = domainwallstate(cmps.CAMPS, rng, N, μ)

function domainwallstate(::Type{T}, rng, N, μ) where {T}
  p = exp(μ)/(exp(μ) + exp(-μ))

  cliff = one(CliffordOperator, N)
  onebitinds = Int[]

  for i in 1:N÷2
    if rand(rng) > p
      push!(onebitinds, i)
    end
  end
  for i in N÷2+1:N
    if rand(rng) < p
      push!(onebitinds, i)
    end
  end

  ψ = computationalcamps(T, N, onebitinds)
  return ψ, onebitinds
end

""" ```computationalcamps([T,] N, onebitinds)```

The N-qubit computational basis state with the qubits in `onebitinds` in |1⟩, as a
CAMPS of type T (`cmps.CAMPS`, the default while the port is staged, or
`MagicTensors.QubitCAMPS`): the MPS is |0…0⟩ and the frame is C = ∏ X_i.
"""
computationalcamps(N, onebitinds) = computationalcamps(cmps.CAMPS, N, onebitinds)

function computationalcamps(::Type{mt.QubitCAMPS}, N, onebitinds)
  C = one(CliffordOperator, N)
  for i in onebitinds
    C[i+N] = -C[i+N] # Z_i ↦ −Z_i
  end
  mps = MPS(ComplexF64, siteinds("Qubit", N), "0")
  return mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(C), false)
end

function computationalcamps(::Type{cmps.CAMPS}, N, onebitinds)
  cliff = one(CliffordOperator, N)
  for i in onebitinds
    cliff[i+N] = -cliff[i+N]
  end

  ψ = cmps.CAMPS(N)
  ψ.Cdag *= inv(cliff)
  return ψ
end

# == Observables ========================================================================

""" ```transferredmagnetization(N, init_ones)```

Return a pp.PauliSum representing the transferred magnetization of Rosenberg et al.,
relative to the initial bitstring with 1s at `init_ones`: ΔN_L − ΔN_R, the change in the
number of 1s in the left half (qubits 1:N÷2) minus that in the right half. It is 0 on the
initial state. When the number of 1s is conserved it equals 2ΔN_L, i.e. it counts +2 for
each 1 that crossed from the right half to the left (−2 for left to right); a single
non-conserving flip changes it by ±1. The factor 2 is Rosenberg et al.'s definition.
"""
function transferredmagnetization(N, init_ones)
  init_zeros = setdiff(1:N, init_ones)

  init_zeros_l = Int[i for i in init_zeros if i <= N÷2]
  init_ones_l = Int[i for i in init_ones if i <= N÷2]
  init_zeros_r = Int[i for i in init_zeros if i > N÷2]
  init_ones_r = Int[i for i in init_ones if i > N÷2]

  Zs_l = [pp.PauliString(N, [:Z], [i], -.5) for i in 1:N÷2]
  Is_l = [pp.PauliString(N, [:I], [i], .5) for i in 1:N÷2]

  Zs_r = [pp.PauliString(N, [:Z], [i], .5) for i in N÷2+1:N]
  Is_r = [pp.PauliString(N, [:I], [i], .5) for i in N÷2+1:N]

  append!(Is_l, [pp.PauliString(N, [:I], [i], -1) for i in init_ones_l])
  append!(Is_r, [pp.PauliString(N, [:I], [i], -1) for i in init_zeros_r])

  obs = pp.PauliSum(N)
  for Zi in Zs_l obs += Zi end
  for Ii in Is_l obs += Ii end
  for Zi in Zs_r obs += Zi end
  for Ii in Is_r obs += Ii end
  return obs
end

# == Circuits ===========================================================================

"""layerends(N, t, circbuilder)

Return a vector of the positions at which each layer ends \
in the circuit specified by the arguments"""
function layerends(N, t, circbuilder)
  layer_depth = length(circbuilder(0,0,1,N)[1])
  layer_ends = collect(layer_depth:layer_depth:t*layer_depth)
  return layer_ends
end