using CampsPP
using Random: seed!

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

# seed!(42)

N = 50
t = 200
χ_campspp = 128
thl_campspp = 1e-15
Nmax_campspp = 200
χ_camps = 256
thl_pp = thl_campspp
Nmax_pp = 100000000

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

output = "output/CircuitDynamics.txt"

evs = campspp_circuit_dynamics(ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; showprogress = true, obsname = "magnetization")

open(output, "a") do f
  println(f, "\n")
end
_ = camps_circuit_dynamics(ψ, 2*χ_camps, gates, phases, obs, output; showprogress = true)

angles = -2 .* phases
open(output, "a") do f
  println(f, "\n")
end
ev = cmps.expectation(ψ, obs)
CampsPP.append_datapoint(output, 0, real(ev))
_ = pauliprop_circuit_dynamics(ψ, 0, thl_pp, gates, angles, Nmax_pp, obs, output; showprogress = true)
return