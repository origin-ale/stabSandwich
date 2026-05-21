using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

# == States =============================================================================

""" ```domainwallstate(N, μ)```

Return an N-qubit product state from the Rosenberg et al. domain wall \
ensemble with contrast parameter μ.
"""

function domainwallstate(N, μ)
  p = exp(μ)/(exp(μ) + exp(-μ))

  states = []
  onebitinds = Int[]

  for i in 1:N÷2
    if rand() < p
      push!(states, "0")
    else
      push!(states, "1")
      push!(onebitinds, i)
    end
  end
  for i in N÷2+1:N
    if rand() < p
      push!(states, "1")
      push!(onebitinds, i)
    else
      push!(states, "0")
    end
  end

  sites = siteinds("Qubit", N)
  ψ = MPS(sites, states)
  ψ = cmps.CAMPS(ψ)

  return ψ, onebitinds
end

# == Observables ========================================================================

""" ```transferredmagnetization(N)```

Return a pp.PauliSum representing the transferred magnetization (Rosenberg et al.) \
ie. the number of 1s which appeared in the left half of the bitstring
"""
function transferredmagnetization(N, onebitinds)
  init_zeros = setdiff(1:N÷2, onebitinds)
  Zs = [pp.PauliString(N, [:Z], [i], -1) for i in init_zeros]
  Is = [pp.PauliString(N, [:I], [i], 1) for i in init_zeros]
  obs = pp.PauliSum(N)
  for Zi in Zs obs += Zi end
  for Ii in Is obs += Ii end
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