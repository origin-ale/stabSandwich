using CampsPP
using Random: seed!

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

# seed!(42)

N = 25
t = 120
χ_campspp = 64
thl_campspp = 1e-13
Nmax_campspp = 200
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = Nmax_campspp

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

output = "output/CircuitDynamics.txt"

printstyled("Running magnetization circuit dynamics until failure for N=$N, χ_campspp = $χ_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)

CampsPP.initialize_output(output, "magnetization", Dict(:N => N, :χ => χ_campspp, :Nmax => Nmax_campspp))
evs = campspp_circuit_dynamics(ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; showprogress = true)

open(output, "a") do f
  println(f, "\n")
end
_ = camps_circuit_dynamics(ψ, gates, phases, χ_camps, obs, output; showprogress = true)

angles = -2 .* phases
open(output, "a") do f
  println(f, "\n")
end
ev = cmps.expectation(ψ, obs)
CampsPP.append_datapoint(output, 0, real(ev))
_ = pauliprop_circuit_dynamics(ψ, 0, gates, angles, thl_pp, Nmax_pp, obs, output; showprogress = true)
return