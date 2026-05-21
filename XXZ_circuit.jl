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
magic_syms = [:X]
magic_inds = [N÷2]
output = "output/XXZDynamics.txt"

layer_depth = length(xxz_circuit(0,0,1,N)[1])
layer_ends = collect(layer_depth:layer_depth:t*layer_depth)

ψ = cmps.CAMPS(N)
onebitinds = Integer[]
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = xxz_circuit(ϕ, θ, t, N)
gates, phases, layer_ends = dopeMagic(N, gates, phases, layer_ends, magic_syms, magic_inds, magic_prob)

printstyled("Running magic-doped XXZ circuit dynamics until failure for N=$N, χ_campspp = $χ_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)

CampsPP.initialize_output(
  output, 
  "Magnetization", 
  Dict(
    "N" => N, 
    "χ" => χ_campspp, 
    "Nmax" => Nmax_campspp))

evs = campspp_circuit_dynamics(
  ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; 
  showprogress = true, layer_ends = layer_ends)

_ = camps_circuit_dynamics(
  ψ, gates, phases, χ_camps, obs, output; 
  showprogress = true, layer_ends = layer_ends)

_ = pauliprop_circuit_dynamics(
  onebitinds, gates, phases, thl_pp, Nmax_pp, obs, output; 
  showprogress = true, layer_ends = layer_ends)

return