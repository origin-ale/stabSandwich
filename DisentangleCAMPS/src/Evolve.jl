using CliffordMPS
using QuantumClifford
using ITensors, ITensorMPS
using ProgressMeter

using Revise

bonddim(ψ::CAMPS) = maximum(dim.(linkinds(ψ.mps)))

generate_showvalues(χ, bd) = () -> [("Bond dimension (max $χ)", bd)]

"```evolve_bonddim(ψ, χ, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along the Pauli rotation circuit\
with Pauli strings paulis and the given phases, until end or bond dim = χ.
Return the evolved CAMPS and stopping time."
function evolve_bonddim(ψ::CAMPS,
                χ::Integer,
                paulistrings::Vector{<:PauliOperator},
                phases::Vector{<:Real};
                showprogress = false)
  k = 0
  s = 0
  progressthresh = ProgressUnknown(0; dt = 0.05, desc = "Evolving CAMPS… t =", enabled = showprogress)
  while s < length(paulistrings) && bonddim(ψ) < χ
    s += 1
    k = apply!(ψ, k, paulistrings[s], phases[s])
    next!(progressthresh; showvalues = generate_showvalues(χ, bonddim(ψ)))
  end
  finish!(progressthresh)
  return ψ, k, s
end

"```evolve(ψ, t, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along the Pauli rotation circuit\
with Pauli strings paulis and the given phases, until layer t.
Return the evolved CAMPS."
function evolve(ψ::CAMPS,
                t::Integer,
                paulistrings::Vector{<:PauliOperator},
                phases::Vector{<:Real};
                showprogress = false)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    k = apply!(ψ, k, paulistrings[s], phases[s])
    next!(progressbar)
  end
  return ψ, k
end

"```evolve_deepcliffords(ψ, t, paulis, phases; [showprogress::Bool])```

Evolve the CAMPS ψ along a Pauli rotation-doped deep Clifford circuit, until layer t.\
Layer s has a 2N²-deep random Clifford and a paulis[s]-rotation with angle phases[s].
Return the evolved CAMPS."
function evolve_deepcliffords(ψ::CAMPS,
                              t::Integer,
                              paulistrings::Vector{<:PauliOperator},
                              phases::Vector{<:Real};
                              showprogress = false)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    applyQOp!(ψ, QOp(:randCliffCircuit))
    k = apply!(ψ, k, paulistrings[s], phases[s])
    next!(progressbar)
  end
  return ψ, k
end

"```apply!(ψ, k, P, ϕ)```

Apply the Pauli operator P to the CAMPS ψ with k free qubits, disentangling if possible.
Modify ψ in-place and return the new number of free qubits."
function CliffordMPS.apply!(ψ::CAMPS,
                            k::Integer,
                            P::PauliOperator,
                            ϕ::Real)
  N = length(ψ)
  I = PauliOperator(0x0, fill(false,N), fill(false,N))
  C = inv(ψ.Cdag)

  aϕ = abs(ϕ)
  if (isapprox(aϕ, 0; atol = 1e-6) || aϕ ≈ π)
    return k
  elseif aϕ ≈ π/4
    apply_pi4_cliff!(ψ, P)
    return k
  elseif aϕ ≈ 3π/4
    apply_3pi4_cliff!(ψ, P)
    return k
  elseif aϕ ≈ π/2
    ψ.Cdag = inv(CliffordOperator(P)) * ψ.Cdag
    return k
  end

  nature = paulinature(k, C, P)
  if nature == :disentanglable
    D, sign = disentangler(k, C, P)
    ψ.Cdag = inv(D) * ψ.Cdag
    addmagicstate!(ψ, k, sign*ϕ)
    k += 1
  elseif nature == :logical
    R = PauliSum([cos(ϕ), im * sin(ϕ)], Stabilizer([I,P]))
    applyGate!(ψ, R)
  elseif nature == :trivial
  end
  return k
end

function apply_pi4_cliff!(ψ::CAMPS, P)
  N = length(ψ)
  op = one(CliffordOperator, N)
  for i in 1:2N
    Q = one(CliffordOperator, N)[i]
    if Bool(comm(P, Q)) # comm(P,Q) == 0x1 if P and Q anticommute
      op[i] = im * P * Q
    end
  end
  ψ.Cdag = ψ.Cdag * inv(op)
end

apply_3pi4_cliff!(ψ::CAMPS, P) = apply_pi4_cliff!(ψ, P)

"```addmagicstate!(ψ, k, phase)```

Turn ψ's (k+1)th qubit from |0⟩ to the Liu and Clark (2025) magic state with given phase"
function addmagicstate!(ψ::CAMPS, k::Integer, phase::Real)
  magifier_os = OpSum()
  magifier_os += cos(phase), "Id", k+1
  magifier_os += im * sin(phase), "X", k+1
  sites = siteinds(ψ.mps)
  magifier = MPO(magifier_os, sites)
  ψ.mps = ITensors.apply(magifier, ψ.mps)
  return nothing
end