using CampsPP
using Random: seed!

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

# seed!(42)

N = 13
t = 30
χ = 64
thl = 1e-10
Nmax = 500

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

evs = campspp_circuit_dynamics(ψ, χ, thl, Nmax, gates, phases, obs, output; showprogress = true, obsname = "magnetization")

ψ = cmps.CAMPS(N)
open(output, "a") do f
  println(f, "\n")
end
_ = camps_circuit_dynamics(ψ, χ, gates, phases, obs, output; showprogress = true)

ψ = cmps.CAMPS(N)
open(output, "a") do f
  println(f, "\n")
end
_ = pauliprop_circuit_dynamics(ψ, 0, thl, gates, phases, Nmax, obs, output; showprogress = true)
return