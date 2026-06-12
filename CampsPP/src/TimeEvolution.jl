using Revise

using DisentangleCAMPS
import PauliPropagation as pp
import CliffordMPS as cmps

using ProgressMeter

# == File manipulation ==================================================================

function initialize_output(output::AbstractString, obsname::AbstractString, params::AbstractDict)
  open(output, "w") do f
    pairs = ["$(string(k))=$(string(v))" for (k,v) in params]
    header = isempty(pairs) ? "" : (join(pairs, " ") * " ")
    println(f, "# $(header) \n# Cycle $obsname")
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

# == Utility ============================================================================

function getlayer(i, layer_ends)
  isnothing(layer_ends) && return i
  for l in eachindex(layer_ends)
    layer_ends[l] >= i && return l
  end
end

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
layer_ends = nothing)
  ψ = copy(ψ_ext)

  ψ_evo, _, s, evs_camps = camps_circuit_dynamics(
    ψ, gates, phases, χ, obs, output; 
    showprogress = showprogress, k = k, layer_ends = layer_ends)

  s, evs_pp = pauliprop_circuit_dynamics(
    ψ_evo, s, gates, phases, thl, Nmax, obs, output; 
    showprogress=showprogress, layer_ends = layer_ends)

  evs_tot = []
  append!(evs_tot, evs_camps)
  append!(evs_tot, evs_pp)

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
layer_ends = nothing)
  return campspp_circuit_dynamics(
    ψ_ext, χ, thl, Nmax, gates, phases, obs, nothing;
    showprogress = showprogress, k = k, layer_ends = layer_ends)
end

function camps_circuit_dynamics(
ψ_ext::cmps.CAMPS,
gates::Vector,
phases::Vector,
χ::Integer,
obs::pp.PauliSum;
showprogress = false,
k = 0,
layer_ends = nothing)
  return camps_circuit_dynamics(
    ψ_ext, gates, phases, χ, obs, nothing;
    showprogress = showprogress, k = k, layer_ends = layer_ends)
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
layer_ends = nothing)

  ψ = copy(ψ_ext)                              
  N = length(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs_camps = []
  progress = ProgressUnknown(desc = "Evolving CAMPS… gate ", enabled = showprogress)
  obs_cmps = cmps.PauliSum(obs)

  append_expectation!(evs_camps, output, ψ, obs_cmps, 0)
  layer += 1

  while DisentangleCAMPS.bonddim(ψ) < χ && i < M
    i += 1
    gate = PauliOperator(getpauli(gates[i], N))
    phase = phases[i]

    k = apply!(ψ, k, gate, phase)

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs_camps, output, ψ, obs_cmps, layer)
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

  return ψ, k, i, evs_camps
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
layer_ends = nothing)
  return pauliprop_circuit_dynamics(
    ψ_ext, start_gate, gates, phases, thl, Nmax, obs, nothing;
    showprogress = showprogress, layer_ends = layer_ends)
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
layer_ends = nothing)

  ψ = copy(ψ_ext)
  M = length(gates)
  i = start_gate
  NP = 1
  layer = getlayer(start_gate+1, layer_ends)
  gates_pp = []
  angles_pp = []
  evs_pp = []
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

  while NP < Nmax && i < M
    i += 1
    gate = gates[i]
    angle = -2. * phases[i]

    push!(gates_pp, gate)
    push!(angles_pp, angle)

    if isnothing(layer_ends) || i == layer_ends[layer]
      paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
      NP = length(paulisum)
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
layer_ends = nothing)
  return pauliprop_circuit_dynamics(
    onebitinds, gates, phases, thl, Nmax, obs, nothing;
    showprogress = showprogress, layer_ends = layer_ends)
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
layer_ends = nothing)

  M = length(gates)
  i = 0
  NP = 1
  layer = 0
  gates_pp = []
  angles_pp = []
  evs_pp = []
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

  append_expectation!(evs_pp, output, onebitinds, obs, 0)
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
  return i, evs_pp
end
