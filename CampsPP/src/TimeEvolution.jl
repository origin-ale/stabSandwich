using Revise
using ProgressMeter

using CampsPP
using DisentangleCAMPS
import PauliPropagation as pp
import CliffordMPS as cmps

function append_datapoint(filename::AbstractString, x::Number, y::Number)
    open(filename, "a") do f
        println(f, "$x\t$y")
    end
end

showvalues_χ(χ, bd) = () -> [("Bond dimension (max $χ)", bd)]
showvalues_P(Pmax, P) = () -> [("Pauli terms (max $Pmax)", P)]

"```camps_rndrotation_dynamics(ψ, χ, obs, output; [showprogress], [k=0])```

Evolve the CAMPS ψ (with k initial magic qubits)
through a random Pauli rotation circuit until bond dim = χ.

At each step, append the expectation value of obs to the output file.

Return the evolved CAMPS and stopping time."
function camps_rndrotation_dynamics(ψ::cmps.CAMPS,
                            χ::Integer,
                            obs::pp.PauliSum,
                            output::AbstractString;
                            showprogress = false,
                            k = 0)
  N = length(ψ)                          
  s = 0
  ev = cmps.expectation(ψ, obs)
  append_datapoint(output, s, real(ev))
  progress = ProgressUnknown(0; desc = "Evolving CAMPS… t =", enabled = showprogress)
  while DisentangleCAMPS.bonddim(ψ) < χ
    s += 1
    gate, phase = random_rotation(N, PauliOperator([1]))
    k = apply!(ψ, k, gate, phase)
    ev = cmps.expectation(ψ, obs)
    append_datapoint(output, s, real(ev))
    next!(progress; showvalues = showvalues_χ(χ, DisentangleCAMPS.bonddim(ψ)))
  end
  finish!(progress)
  return ψ, k, s
end

function pauliprop_rndrotation_dynamics(ψ::cmps.CAMPS,
                                        s0::Integer,
                                        thl::Real,
                                        Nmax::Integer,
                                        obs::pp.PauliSum,
                                        output::AbstractString;
                                        showprogress = false)
  N = length(ψ)
  s = s0
  NP = 1
  gates_pp = []
  angles_pp = []
  progress = ProgressUnknown(dt = 0.05, desc = "Evolving with Pauli prop… t =", enabled = showprogress)
  while NP < Nmax
    s += 1
    gate, angle = random_rotation(N, pp.PauliRotation([:X],[1]))
    push!(gates_pp, gate)
    push!(angles_pp, angle)
    paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
    NP = length(paulisum)
    ev = cmps.expectation(ψ, paulisum)
    append_datapoint(output, s, real(ev))
    update!(progress, s; showvalues = showvalues_P(Nmax, NP))
  end
  finish!(progress)
  return s
end