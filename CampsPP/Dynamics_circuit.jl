using CampsPP
using Random: seed!

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf

seed!(42)

N = 13
t = 30
χ = 64
thl = 1e-10
Nmax = 200

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

ψ_evo, _, s, evs_camps = campspp_circuit_dynamics(ψ, χ, thl, Nmax, gates, phases, obs, output; showprogress = true, obsname = "magnetization")