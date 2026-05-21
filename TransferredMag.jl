using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!

seed!(2)

N = 46
t = 23
ϕ = π/4
θ = π/4

χ_campspp = 64
thl_campspp = 1e-15
Nmax_campspp = 200
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = Nmax_campspp
magic_prob = 0
magic_syms = [:X, :X]
magic_inds = [N÷2, N÷2+1]
output = "output/TMDynamics.txt"

layer_depth = length(xxz_circuit(0,0,1,N)[1])
layer_ends = collect(layer_depth:layer_depth:t*layer_depth)

sites = siteinds("Qubit", N)
states = [i<=N÷2 ? "0" : "1" for i in eachindex(sites)]
ψ = MPS(sites, states)
ψ = cmps.CAMPS(ψ)
onebitinds = collect(N÷2+1:N)

Zs = [pp.PauliString(N, [:Z], [i], -1) for i = 1:N÷2]
Is = [pp.PauliString(N, [:I], [i], 1) for i = 1:N÷2]
obs = pp.PauliSum(N)
for Zi in Zs global obs += Zi end
for Ii in Is global obs += Ii end

gates, phases = xxz_circuit(ϕ, θ, t, N)
gates, phases, layer_ends = dopeMagic(N, gates, phases, layer_ends, magic_syms, magic_inds, magic_prob)

printstyled("Running magic-doped XXZ circuit dynamics until t = $t for N=$N, χ_campspp = $χ_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)

initialize_output(
  output, 
  "Transferred magnetization", 
  Dict(
    "N" => N, 
    "χ" => χ_campspp, 
    "Nmax" => Nmax_campspp))

evs = campspp_circuit_dynamics(
  ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; 
  showprogress = true, layer_ends = layer_ends)

evs = pauliprop_circuit_dynamics(
  onebitinds, gates, phases, thl_campspp, Nmax_campspp, obs, output; 
  showprogress = true, layer_ends = layer_ends)

return