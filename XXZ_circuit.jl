using Revise
using CampsPP
using DisentangleCAMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!

seed!(2)

N = 22
t = 10
ϕ = π/4
θ = π/4

χ_campspp = 64
thl_campspp = 1e-15
Nmax_campspp = 200
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = Nmax_campspp
magic_prob = .5
end_xxz = [:Z, :Z]
output = "output/XXZDynamics.txt"

ψ = cmps.CAMPS(N)
onebitinds = Integer[]
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = xxz_circuit(ϕ, θ, t, N)
gates, phases, magic_pos = dopeT(N, gates, phases, magic_prob)

printstyled("Running magic-doped XXZ circuit dynamics until failure for N=$N, χ_campspp = $χ_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)
println("Magic doping at $(length(magic_pos)) gates of $(length(gates)).")

CampsPP.initialize_output(
  output, 
  "Magnetization", 
  Dict(
    "N" => N, 
    "χ" => χ_campspp, 
    "Nmax" => Nmax_campspp))

evs = campspp_circuit_dynamics(
  ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; 
  showprogress = true, ev_at = end_xxz)

_ = camps_circuit_dynamics(
  ψ, gates, phases, χ_camps, obs, output; 
  showprogress = true, ev_at = end_xxz)

_ = pauliprop_circuit_dynamics(
  onebitinds, gates, phases, thl_pp, Nmax_pp, obs, output; 
  showprogress = true, ev_at = end_xxz)

return