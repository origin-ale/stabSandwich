using Revise
using ProgressMeter

using CampsPP
using DisentangleCAMPS
import PauliPropagation as pp
import CliffordMPS as cmps

# -- Utility ------------------------------------------------------------------

function append_datapoint(filename::AbstractString, x::Number, y::Number)
    open(filename, "a") do f
        println(f, "$x\t$y")
    end
end

showvalues_χ(χ, bd, N, k) = () -> [("Bond dimension (max $χ)", bd), ("Magic qubits (max $N)", k)]
showvalues_P(Pmax, P) = () -> [("Pauli terms (max $Pmax)", P)]

# -- Random evolution ---------------------------------------------------------

"```camps_rndrotation_dynamics(ψ, χ, obs, output; [showprogress], [k=0])```

Evolve the CAMPS ψ (with k initial magic qubits)
through a random Pauli rotation circuit until bond dim = χ.

At each step, append the expectation value of obs to the output file.

Return the evolved CAMPS, stopping time and vector of expectation values."
function camps_rndrotation_dynamics(ψ_ext::cmps.CAMPS,
                            χ::Integer,
                            obs::pp.PauliSum,
                            output::AbstractString;
                            showprogress = false,
                            k = 0)
  ψ = copy(ψ_ext)
  N = length(ψ)                          
  s = 0
  evs_camps = []
  ev = cmps.expectation(ψ, obs)
  push!(evs_camps, real(ev))
  append_datapoint(output, s, real(ev))
  progress = ProgressUnknown(0; desc = "Evolving CAMPS… gate", enabled = showprogress)
  while DisentangleCAMPS.bonddim(ψ) < χ
    s += 1
    gate, phase = random_rotation(N, PauliOperator([1]))
    k = apply!(ψ, k, gate, phase)
    ev = cmps.expectation(ψ, obs)
    push!(evs_camps, real(ev))
    append_datapoint(output, s, real(ev))
    next!(progress; showvalues = showvalues_χ(χ, DisentangleCAMPS.bonddim(ψ), N, k))
  end
  finish!(progress)
  return ψ, k, s, evs_camps
end

"```pauliprop_rndrotation_dynamics(ψ, s0, thl, Nmax, obs, output; [showprogress])```

Track evolution of the expectation value of obs starting at time s0,
through a random Pauli rotation circuit until its Heisenberg evolution 
contains Nmax Pauli terms.

Use Pauli propagation with coefficient truncation at thl.

At each step, append the expectation value to the output file.

Return the stopping time and vector of expectation values."
function pauliprop_rndrotation_dynamics(ψ_ext::cmps.CAMPS,
                                        s0::Integer,
                                        thl::Real,
                                        Nmax::Integer,
                                        obs::pp.PauliSum,
                                        output::AbstractString;
                                        showprogress = false)
  ψ = copy(ψ_ext)
  N = length(ψ)
  s = s0
  NP = 1
  evs_pp = []
  gates_pp = []
  angles_pp = []
  progress = ProgressUnknown(dt = 0.05, desc = "Evolving with Pauli prop… gate", enabled = showprogress)
  while NP < Nmax
    s += 1
    gate, angle = random_rotation(N, pp.PauliRotation([:X],[1]))
    push!(gates_pp, gate)
    push!(angles_pp, angle)
    paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
    NP = length(paulisum)
    ev = cmps.expectation(ψ, paulisum)
    push!(evs_pp, real(ev))
    append_datapoint(output, s, real(ev))
    update!(progress, s; showvalues = showvalues_P(Nmax, NP))
  end
  finish!(progress)
  return s, evs_pp
end

# -- Circuit evolution --------------------------------------------------------

function camps_circuit_dynamics(ψ_ext::cmps.CAMPS,
                                χ::Integer,
                                gates::Vector,
                                phases::Vector,
                                obs::pp.PauliSum,
                                output::AbstractString;
                                showprogress = false,
                                k = 0,
                                ev_at = nothing)
  ψ = copy(ψ_ext)                              
  N = length(ψ)
  T = length(gates)
  s = 0
  evs_camps = []
  ev = cmps.expectation(ψ, obs)
  push!(evs_camps, real(ev))
  append_datapoint(output, s, real(ev))
  progress = ProgressUnknown(0; desc = "Evolving CAMPS… gate ", enabled = showprogress)
  obs_cmps = cmps.PauliSum(obs)
  while DisentangleCAMPS.bonddim(ψ) < χ && s < T
    s += 1
    gate = PauliOperator(getpauli(gates[s], N))
    phase = phases[s]
    k = apply!(ψ, k, gate, phase)
    if isnothing(ev_at) || gates[s].symbols == ev_at
      ev = cmps.expectation(ψ, obs_cmps)
      push!(evs_camps, real(ev))
      append_datapoint(output, s, real(ev))
    end
    next!(progress; showvalues = showvalues_χ(χ, DisentangleCAMPS.bonddim(ψ), N, k))
  end
  finish!(progress)
  return ψ, k, s, evs_camps
end

function pauliprop_circuit_dynamics(ψ_ext::cmps.CAMPS,
                                    s0::Integer,
                                    thl::Real,
                                    gates::Vector,
                                    angles::Vector,
                                    Nmax::Integer,
                                    obs::pp.PauliSum,
                                    output::AbstractString;
                                    showprogress = false,
                                    ev_at = nothing)
  ψ = copy(ψ_ext)
  Tp = length(gates)
  s = s0
  NP = 1
  evs_pp = []
  gates_pp = []
  angles_pp = []
  progress = ProgressUnknown(dt = 0.05, desc = "Evolving with Pauli prop… gate ", enabled = showprogress)
  while NP < Nmax && s-s0 < Tp
    s += 1
    gate = gates[s-s0]
    angle = angles[s-s0]
    push!(gates_pp, gate)
    push!(angles_pp, angle)
    paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
    NP = length(paulisum)
    if isnothing(ev_at) || gates[s].symbols == ev_at
      if s0 == 0
        ev = pp.overlapwithzero(paulisum)
      else
        ev = cmps.expectation(ψ, paulisum)
      end
      push!(evs_pp, real(ev))
      append_datapoint(output, s, real(ev))
    end
    update!(progress, s; showvalues = showvalues_P(Nmax, NP))
  end
  finish!(progress)
  return s, evs_pp
end

function campspp_circuit_dynamics(ψ_ext::cmps.CAMPS,
                                  χ::Integer,
                                  thl::Real,
                                  Nmax::Integer,
                                  gates::Vector{},
                                  phases::Vector,
                                  obs::pp.PauliSum,
                                  output::AbstractString;
                                  showprogress = false,
                                  k = 0,
                                  obsname::AbstractString = "[unknown]",
                                  ev_at = nothing)
  ψ = copy(ψ_ext)
  N = length(ψ)
  open(output, "w") do f
    println(f, "# N=$N χ=$χ Nmax=$Nmax obs=$obsname")
  end

  ψ_evo, _, s, evs_camps = camps_circuit_dynamics(ψ, χ, gates, phases, obs, output; showprogress = showprogress, k = k, ev_at = ev_at)
  open(output, "a") do f
    println(f, "# CAMPS stopped at gate $s (χ ≥ $χ)")
  end

  leftover_gates, leftover_angles = leftover_rotgates(s, gates, phases)
  s, evs_pp = pauliprop_circuit_dynamics(ψ_evo, s, thl, leftover_gates, leftover_angles, Nmax, obs, output; showprogress=showprogress, ev_at = ev_at)

  open(output, "a") do f
    println(f, "# Pauli prop. stopped at gate $s (N_pauli ≥ $Nmax)")
  end

  evs_tot = []
  append!(evs_tot, evs_camps)
  append!(evs_tot, evs_pp)
  return evs_tot
end