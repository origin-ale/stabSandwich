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
out_array::Vector, out_file::AbstractString, onebitinds::Vector{<:Integer}, obs, i::Integer)
  ev = pp.overlapwithcomputational(obs, onebitinds)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, ev::Number, i::Integer)
  push!(out_array, real(ev))
  append_datapoint(out_file, i, real(ev))
end

# == ProgressMeter showvalues ===========================================================

showvalues_χ(χ, bd, N, k) = () -> [("Bond dimension (max $χ)", bd), ("Magic qubits (max $N)", k)]
showvalues_P(Pmax, P) = () -> [("Pauli terms (max $Pmax)", P)]

# == Circuit evolution ==================================================================

function campspp_circuit_dynamics(
ψ_ext::cmps.CAMPS,
χ::Integer,
thl::Real,
Nmax::Integer,
gates::Vector{},
phases::Vector,
obs::pp.PauliSum,
output::AbstractString;
showprogress = false,
k = 0,
ev_at = nothing)
  ψ = copy(ψ_ext)

  ψ_evo, _, s, evs_camps = camps_circuit_dynamics(
    ψ, gates, phases, χ, obs, output; 
    showprogress = showprogress, k = k, ev_at = ev_at)

  s, evs_pp = pauliprop_circuit_dynamics(
    ψ_evo, s, gates, phases, thl, Nmax, obs, output; 
    showprogress=showprogress, ev_at = ev_at)

  evs_tot = []
  append!(evs_tot, evs_camps)
  append!(evs_tot, evs_pp)

  return evs_tot
end

function camps_circuit_dynamics(
ψ_ext::cmps.CAMPS,
gates::Vector,
phases::Vector,
χ::Integer,
obs::pp.PauliSum,
output::AbstractString;
showprogress = false,
k = 0,
ev_at = nothing)

  ψ = copy(ψ_ext)                              
  N = length(ψ)
  M = length(gates)
  i = 0
  evs_camps = []
  progress = ProgressUnknown(desc = "Evolving CAMPS… gate ", enabled = showprogress)
  obs_cmps = cmps.PauliSum(obs)

  append_expectation!(evs_camps, output, ψ, obs_cmps, 0)

  while DisentangleCAMPS.bonddim(ψ) < χ && i < M
    i += 1
    gate = PauliOperator(getpauli(gates[i], N))
    phase = phases[i]

    k = apply!(ψ, k, gate, phase)

    if isnothing(ev_at) || gates[i].symbols == ev_at
      append_expectation!(evs_camps, output, ψ, obs_cmps, i)
    end

    bd = DisentangleCAMPS.bonddim(ψ)
    next!(progress; showvalues = showvalues_χ(χ, bd, N, k))
  end

  finish!(progress)
  open(output, "a") do f
    reason = (i == M) ? "(end of circuit)" : "(bond dim. ≥ $χ)"
    println(f, "# CAMPS stopped at gate $i ", reason, "\n\n")
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
obs::pp.PauliSum,
output::AbstractString;
showprogress = false,
ev_at = nothing)

  ψ = copy(ψ_ext)
  M = length(gates)
  i = start_gate
  NP = 1
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
    paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)

    NP = length(paulisum)
    if isnothing(ev_at) || gates[i].symbols == ev_at
      append_expectation!(evs_pp, output, ψ, paulisum, i)
    end

    next!(progress; showvalues = showvalues_P(Nmax, NP))
  end

  finish!(progress)
  open(output, "a") do f
    reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
    println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
  end
  return i, evs_pp
end

function pauliprop_circuit_dynamics(
onebitinds::Vector{<:Integer},
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum,
output::AbstractString;
showprogress = false,
ev_at = nothing)

  M = length(gates)
  i = 0
  NP = 1
  gates_pp = []
  angles_pp = []
  evs_pp = []
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

  append_expectation!(evs_pp, output, onebitinds, obs, 0)

  while NP < Nmax && i < M
    i += 1
    gate = gates[i]
    angle = -2. * phases[i]

    push!(gates_pp, gate)
    push!(angles_pp, angle)
    paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)

    NP = length(paulisum)
    if isnothing(ev_at) || gates[i].symbols == ev_at
      append_expectation!(evs_pp, output, onebitinds, paulisum, i)
    end

    next!(progress; showvalues = showvalues_P(Nmax, NP))
  end

  finish!(progress)
  open(output, "a") do f
    reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
    println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
  end
  return i, evs_pp
end