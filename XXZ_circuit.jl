using Revise
using CampsPP
using DisentangleCAMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!

seed!(2)

N = 50
t = 100
χ_campspp = 64
thl_campspp = 1e-15
Nmax_campspp = 200
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = Nmax_campspp
magic_prob = .5

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

magic_pos = []
gates, phases = CampsPP.xxz_circuit(t, N)
phases = [π/4 for g in gates]
newgates = copy(gates)
newphases = copy(phases)
os = 0
for i in eachindex(gates)
  if rand() < magic_prob
    insert!(newgates, i+os, pp.PauliRotation([:Z], rand(1:N)))
    insert!(newphases, i+os, π/8)
    push!(magic_pos, i+os)
    global os += 1
  end
end
gates = newgates
phases = newphases

output = "output/XXZDynamics.txt"

printstyled("Running magic-doped XXZ circuit dynamics until failure for N=$N, χ_campspp = $χ_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)
println("Magic doping at $(length(magic_pos)) gates of $(length(gates)).")

evs = campspp_circuit_dynamics(ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; showprogress = true, obsname = "magnetization")

open(output, "a") do f
  println(f, "\n")
end
_ = camps_circuit_dynamics(ψ, χ_camps, gates, phases, obs, output; showprogress = true)

angles = -2 .* phases
open(output, "a") do f
  println(f, "\n")
end
ev = cmps.expectation(ψ, obs)
CampsPP.append_datapoint(output, 0, real(ev))
_ = pauliprop_circuit_dynamics(ψ, 0, thl_pp, gates, angles, Nmax_pp, obs, output; showprogress = true)
return