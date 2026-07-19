using Revise

using DisentangleCAMPS
import PauliPropagation as pp
import CliffordMPS as cmps
using ITensors, ITensorMPS

using ProgressMeter

# == File manipulation ==================================================================

function initialize_output(output::AbstractString, obsname::AbstractString, params::AbstractDict)
  open(output, "w") do f
    pairs = ["$(string(k))=$(string(v))" for (k,v) in params]
    header = isempty(pairs) ? "" : (join(pairs, " ") * " ")
    println(f, "# $(header) \n# Cycle $obsname\n\n")
  end
end

function append_datapoint(output::AbstractString, x::Number, y::Number)
    open(output, "a") do f
        println(f, "$x\t$y")
    end
end

# == Append expectation values ==========================================================

function append_expectation!(
out_array::Vector, out_file::AbstractString, ψ::cmps.CAMPS, obs, i::Integer)
  ev = cmps.expectation(ψ, obs)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, ::Nothing, ψ::cmps.CAMPS, obs, i::Integer)
  ev = cmps.expectation(ψ, obs)
  push!(out_array, real(ev))
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, ψ::ITensorMPS.MPS, obs::cmps.PauliSum, i::Integer)
  ev = cmps.expectation(ψ, obs)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, ::Nothing, ψ::ITensorMPS.MPS, obs::cmps.PauliSum, i::Integer)
  ev = cmps.expectation(ψ, obs)
  push!(out_array, real(ev))
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, onebitinds::Vector{<:Integer}, obs, i::Integer)
  ev = pp.overlapwithcomputational(obs, onebitinds)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, ::Nothing, onebitinds::Vector{<:Integer}, obs, i::Integer)
  ev = pp.overlapwithcomputational(obs, onebitinds)
  push!(out_array, real(ev))
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, ev::Number, i::Integer)
  push!(out_array, real(ev))
  append_datapoint(out_file, i, real(ev))
end

function append_expectation!(
out_array::Vector, ::Nothing, ev::Number, i::Integer)
  push!(out_array, real(ev))
end

# == ProgressMeter showvalues ===========================================================

showvalues_χ(χ, bd, N, k) = () -> [("Bond dimension (max $χ)", bd), ("Magic qubits (max $N)", k)]
showvalues_P(Pmax, P) = () -> [("Pauli terms (max $Pmax)", P)]
showvalues_ee(ee) = () -> [("Entanglement entropies", round.(ee; digits=2))]

# == Utility ============================================================================

function getlayer(i, layer_ends)
  isnothing(layer_ends) && return i
  for l in eachindex(layer_ends)
    layer_ends[l] >= i && return l
  end
end

meanweight(psum::pp.PauliSum) = mean(pp.countweight(psum))

# == Circuit evolution ==================================================================

function campspp_circuit_dynamics(
ψ_ext::cmps.CAMPS,
χ::Integer,
thl::Real,
Nmax::Integer,
gates::Vector{},
phases::Vector,
obs::pp.PauliSum,
output;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)
  ψ = deepcopy(ψ_ext)

  if track
    ψ_evo, _, s, evs_camps, bonddims = camps_circuit_dynamics(
      ψ, gates, phases, χ, obs, output;
      showprogress = showprogress, k = k, layer_ends = layer_ends, track = true)

    s, evs_pp, nterms, avgweights = pauliprop_circuit_dynamics(
      ψ_evo, s, gates, phases, thl, Nmax, obs, output;
      showprogress = showprogress, layer_ends = layer_ends, track = true)
  else
    ψ_evo, _, s, evs_camps = camps_circuit_dynamics(
      ψ, gates, phases, χ, obs, output;
      showprogress = showprogress, k = k, layer_ends = layer_ends)

    s, evs_pp = pauliprop_circuit_dynamics(
      ψ_evo, s, gates, phases, thl, Nmax, obs, output;
      showprogress=showprogress, layer_ends = layer_ends)
  end

  evs_tot = []
  append!(evs_tot, evs_camps)
  append!(evs_tot, evs_pp)

  track && return evs_tot, s, bonddims, nterms, avgweights
  return evs_tot, s
end

function campspp_circuit_dynamics(
ψ_ext::cmps.CAMPS,
χ::Integer,
thl::Real,
Nmax::Integer,
gates::Vector{},
phases::Vector,
obs::pp.PauliSum;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)
  return campspp_circuit_dynamics(
    ψ_ext, χ, thl, Nmax, gates, phases, obs, nothing;
    showprogress = showprogress, k = k, layer_ends = layer_ends, track = track)
end

function camps_circuit_dynamics(
ψ_ext::cmps.CAMPS,
gates::Vector,
phases::Vector,
χ::Integer,
obs::pp.PauliSum;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)
  return camps_circuit_dynamics(
    ψ_ext, gates, phases, χ, obs, nothing;
    showprogress = showprogress, k = k, layer_ends = layer_ends, track = track)
end

function camps_circuit_dynamics(
ψ_ext::cmps.CAMPS,
gates::Vector,
phases::Vector,
χ::Integer,
obs::pp.PauliSum,
output;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)

  ψ = deepcopy(ψ_ext)
  N = length(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs_camps = []
  bonddims = Int[]
  progress = ProgressUnknown(desc = "Evolving CAMPS… gate ", enabled = showprogress)
  obs_cmps = cmps.PauliSum(obs)

  append_expectation!(evs_camps, output, ψ, obs_cmps, 0)
  track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
  layer += 1

  while DisentangleCAMPS.bonddim(ψ) < χ && i < M
    i += 1
    gate = PauliOperator(getpauli(gates[i], N))
    phase = phases[i]

    k = apply!(ψ, k, gate, phase)

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs_camps, output, ψ, obs_cmps, layer)
      track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
      layer += 1
    end

    bd = DisentangleCAMPS.bonddim(ψ)
    next!(progress; showvalues = showvalues_χ(χ, bd, N, k))
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(bond dim. ≥ $χ)"
      println(f, "# CAMPS stopped at gate $i ", reason, "\n\n")
    end
  end

  track && return ψ, k, i, evs_camps, bonddims
  return ψ, k, i, evs_camps
end

function campssrc_circuit_dynamics(
ψ_ext::cmps.CAMPS,
gates::Vector,
phases::Vector,
thl::Real,
obs::pp.PauliSum;
criterion = :entangle,
strategy = :brickwork,
showprogress = false,
layer_ends = nothing,
track = false)
  return campssrc_circuit_dynamics(
    ψ_ext, gates, phases, thl, obs, nothing;
    criterion = criterion, strategy = strategy,
    showprogress = showprogress, layer_ends = layer_ends, track = track)
end

function campssrc_circuit_dynamics(
ψ_ext::cmps.CAMPS,
gates::Vector,
phases::Vector,
thl::Real,
obs::pp.PauliSum,
output;
criterion = :entangle,
strategy = :brickwork,
showprogress = false,
layer_ends = nothing,
track = false)

  ψ = deepcopy(ψ_ext)
  N = length(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs_camps = []
  bonddims = Int[]
  progress = ProgressUnknown(desc = "Evolving CAMPS (SVD)… gate ", enabled = showprogress)
  obs_cmps = cmps.PauliSum(obs)
  crit = cmps.DisentangleCriterion(criterion)
  strat = cmps.DisentangleStrategy(strategy)

  # Gates are unitary and disentangling Cliffords are exact, so the norm only
  # shrinks when the SVD cutoff in apply! actually discards weight.
  norm0_sq = norm(ψ)^2
  truncated = false

  append_expectation!(evs_camps, output, ψ, obs_cmps, 0)
  track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
  layer += 1

  while i < M
    i += 1
    gate = PauliOperator(getpauli(gates[i], N))
    phase = phases[i]

    apply!(ψ, gate, phase, thl)

    if !truncated && norm(ψ)^2 < norm0_sq * (1 - thl)
      truncated = true
    end

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs_camps, output, ψ, obs_cmps, layer)
      truncated || cmps.disentangle!(ψ, strat, N; criterion = crit, min_diff = 6)
      track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
      layer += 1
    end

    bd = DisentangleCAMPS.bonddim(ψ)
    next!(progress; showvalues = () -> [("Bond dimension", bd)])
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      println(f, "# CAMPS (SVD) stopped at gate $i (end of circuit)\n\n")
    end
  end

  track && return ψ, i, evs_camps, bonddims
  return ψ, i, evs_camps
end

function pauliprop_circuit_dynamics(
ψ_ext::cmps.CAMPS,
start_gate::Integer,
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum;
showprogress = false,
layer_ends = nothing,
track = false)
  return pauliprop_circuit_dynamics(
    ψ_ext, start_gate, gates, phases, thl, Nmax, obs, nothing;
    showprogress = showprogress, layer_ends = layer_ends, track = track)
end

function pauliprop_circuit_dynamics(
ψ_ext::cmps.CAMPS,
start_gate::Integer,
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum,
output;
showprogress = false,
layer_ends = nothing,
track = false)

  ψ = deepcopy(ψ_ext)
  M = length(gates)
  i = start_gate
  NP = 1
  layer = getlayer(start_gate+1, layer_ends)
  gates_pp = []
  angles_pp = []
  evs_pp = []
  nterms = Int[]
  avgweights = Float64[]
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

  track && push!(nterms, length(obs))
  track && push!(avgweights, meanweight(obs))

  while NP < Nmax && i < M
    i += 1
    gate = gates[i]
    angle = -2. * phases[i]

    push!(gates_pp, gate)
    push!(angles_pp, angle)

    if isnothing(layer_ends) || i == layer_ends[layer]
      paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
      NP = length(paulisum)
      track && push!(nterms, NP)
      track && push!(avgweights, meanweight(paulisum))
      append_expectation!(evs_pp, output, ψ, paulisum, layer)
      layer += 1
    end

    next!(progress; showvalues = showvalues_P(Nmax, NP))
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
      println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
    end
  end
  track && return i, evs_pp, nterms, avgweights
  return i, evs_pp
end

function pauliprop_circuit_dynamics(
onebitinds::Vector{<:Integer},
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum;
showprogress = false,
layer_ends = nothing,
track = false)
  return pauliprop_circuit_dynamics(
    onebitinds, gates, phases, thl, Nmax, obs, nothing;
    showprogress = showprogress, layer_ends = layer_ends, track = track)
end

function pauliprop_circuit_dynamics(
onebitinds::Vector{<:Integer},
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum,
output::Union{AbstractString, Nothing};
showprogress = false,
layer_ends = nothing,
track = false)

  M = length(gates)
  i = 0
  NP = 1
  layer = 0
  gates_pp = []
  angles_pp = []
  evs_pp = []
  nterms = Int[]
  avgweights = Float64[]
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

  append_expectation!(evs_pp, output, onebitinds, obs, 0)
  track && push!(nterms, length(obs))
  track && push!(avgweights, meanweight(obs))
  layer += 1

  while NP < Nmax && i < M
    i += 1
    gate = gates[i]
    angle = -2. * phases[i]

    push!(gates_pp, gate)
    push!(angles_pp, angle)

    if isnothing(layer_ends) || i == layer_ends[layer]
      paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
      NP = length(paulisum)
      track && push!(nterms, NP)
      track && push!(avgweights, meanweight(paulisum))
      append_expectation!(evs_pp, output, onebitinds, paulisum, layer)
      layer += 1
    end

    next!(progress; showvalues = showvalues_P(Nmax, NP))
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
      println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
    end
  end
  track && return i, evs_pp, nterms, avgweights
  return i, evs_pp
end

function mps_circuit_dynamics(
  ψ_ext::MPS,
  gates::Vector,
  phases::Vector,
  thl,
  obs::pp.PauliSum,
  output;
  showprogress = false,
  layer_ends = nothing
)
  ψ = deepcopy(ψ_ext)                              
  N = length(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs = []
  progress = ProgressUnknown(desc = "Evolving MPS… gate ", enabled = showprogress)
  obs_cmps = cmps.PauliSum(obs)

  append_expectation!(evs, output, ψ, obs_cmps, 0)
  layer += 1

  while i < M
    i += 1
    gate_pauli = PauliOperator(getpauli(gates[i], N))
    I = PauliOperator(0x0, fill(false,N), fill(false,N))
    phase = phases[i]
    gate = cmps.PauliSum([cos(phase), im * sin(phase)], Stabilizer([I,gate_pauli]))

    ψ = ITensorMPS.apply(ψ, gate; cutoff = thl)
    normalize!(ψ)

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs, output, ψ, obs_cmps, layer)
      layer += 1
    end
    ee = cmps.eEntropys!(ψ)
    next!(progress; showvalues = showvalues_ee(ee))
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = "end of circuit"
      println(f, "# MPS stopped at gate $i ", reason, "\n\n")
    end
  end

  return ψ, evs
end