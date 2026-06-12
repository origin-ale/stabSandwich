using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!
using Strided
using LinearAlgebra

Strided.disable_threads()
BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

seed!(2)

N = 10
t = 100
θ = π/2
ϕ = π

χ_campspp = 64
thl_campspp = 1e-15
Nmax_campspp = 200
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = Nmax_campspp
magic_prob = 0.5
warn_on_prestop = true

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

layer_ends = layerends(N, t, fSim_circuit)
gates, phases = fSim_circuit(θ, ϕ, t, N)
n_undoped = length(gates)
gates, phases, layer_ends = dopeT(N, gates, phases, layer_ends, magic_prob)

output = "output/fSimDynamics_$(N).txt"
param_info = Dict(
  "N" => N,
  "t" => t,
  "θ" => θ,
  "ϕ" => ϕ,
  "magic_prob" => magic_prob,
  "χ_campspp" => χ_campspp,
  "thl_campspp" => thl_campspp,
  "Nmax_campspp" => Nmax_campspp,
  "χ_camps" => χ_camps,
  "thl_pp" => thl_pp,
  "Nmax_pp" => Nmax_pp,
  "n.gates" => length(phases))

printstyled("Running magic-doped fSim circuit dynamics until failure for N=$N, χ_campspp = $χ_campspp, thl_campspp = $thl_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)
println("Magic doping at $(length(gates) - n_undoped) gates of $(length(gates)).")

initialize_output(output, "magnetization", param_info)
_, t_stop = campspp_circuit_dynamics(
  ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output;
  showprogress = true, layer_ends = layer_ends)
if warn_on_prestop && (t_stop < length(phases))
  printstyled("WARNING: CAMPS+PP run stopped at gate $t_stop/$(length(phases))\n"; color = :yellow)
end

_ = camps_circuit_dynamics(
  ψ, gates, phases, χ_camps, obs, output;
  showprogress = true, layer_ends = layer_ends)

ev = cmps.expectation(ψ, obs)
append_datapoint(output, 0, real(ev))
_ = pauliprop_circuit_dynamics(
  ψ, 0, gates, phases, thl_pp, Nmax_pp, obs, output;
  showprogress = true, layer_ends = layer_ends)

return
