using CampsPP
using Random: seed!

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

# seed!(42)

N = 21
t = 50
χ = 64
thl = 1e-10
Nmax = 200

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

output = "output/CircuitDynamics.txt"

evs = campspp_circuit_dynamics(ψ, χ, thl, Nmax, gates, phases, obs, output; showprogress = true, obsname = "magnetization")

open(output, "a") do f
  println(f, "\n")
end
_ = camps_circuit_dynamics(ψ, 2*χ, gates, phases, obs, output; showprogress = true)

angles = -2 .* phases
open(output, "a") do f
  println(f, "\n")
end
ev = cmps.expectation(ψ, obs)
CampsPP.append_datapoint(output, 0, real(ev))
_ = pauliprop_circuit_dynamics(ψ, 0, thl, gates, angles, 2*Nmax, obs, output; showprogress = true)
return