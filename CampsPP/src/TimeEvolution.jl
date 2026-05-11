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

generate_showvalues(χ, bd) = () -> [("Bond dimension (max $χ)", bd)]

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
  progressthresh = ProgressUnknown(0; dt = 0.05, desc = "Evolving CAMPS… t =", enabled = showprogress)
  while DisentangleCAMPS.bonddim(ψ) < χ
    s += 1
    gate, phase = random_rotation(N)
    k = apply!(ψ, k, gate, phase)
    next!(progressthresh; showvalues = generate_showvalues(χ, DisentangleCAMPS.bonddim(ψ)))
    ev = cmps.expectation(ψ, obs)
    append_datapoint(output, s, real(ev))
  end
  finish!(progressthresh)
  return ψ, k, s
end