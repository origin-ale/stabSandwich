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

  cliff = one(CliffordOperator, N)
  onebitinds = Int[]

  for i in 1:N÷2
    if rand() > p
      push!(onebitinds, i)
    end
  end
  for i in N÷2+1:N
    if rand() < p
      push!(onebitinds, i)
    end
  end

  ψ = computationalcamps(N, onebitinds)
  return ψ, onebitinds
end

function computationalcamps(N, onebitinds)
  cliff = one(CliffordOperator, N)
  for i in onebitinds
    cliff[i+N] = -cliff[i+N]
  end

  ψ = cmps.CAMPS(N)
  ψ.Cdag *= inv(cliff)
  return ψ
end

# == Observables ========================================================================

""" ```transferredmagnetization(N)```

Return a pp.PauliSum representing the transferred magnetization (Rosenberg et al.) \
ie. the number of 1s which appeared in the left half of the bitstring
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