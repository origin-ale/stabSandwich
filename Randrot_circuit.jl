using Revise
using CampsPP
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

# seed!(42)

N = 25
t = 120
χ_campspp = 64
thl_campspp = 1e-13
Nmax_campspp = 200
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = Nmax_campspp
warn_on_prestop = true

ψ = cmps.CAMPS(N)
Zs = [pp.PauliString(N, [:Z], [i], 1/N) for i = 1:N]
obs = pp.PauliSum(Zs)

gates, phases = rotation_circuit(t, N)

output = "output/CircuitDynamics_$(N).txt"
param_info = Dict(
  "N" => N,
  "χ_campspp" => χ_campspp,
  "thl_campspp" => thl_campspp,
  "Nmax_campspp" => Nmax_campspp,
  "χ_camps" => χ_camps,
  "thl_pp" => thl_pp,
  "Nmax_pp" => Nmax_pp,
  "n.gates" => length(phases))

printstyled("Running magnetization circuit dynamics until failure for N=$N, χ_campspp = $χ_campspp, thl_campspp = $thl_campspp, Nmax_campspp = $Nmax_campspp.\n"; color = :cyan)

initialize_output(output, "magnetization", param_info)
_, t_stop = campspp_circuit_dynamics(ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output; showprogress = true)
if warn_on_prestop && (t_stop < length(phases))
  printstyled("WARNING: CAMPS+PP run stopped at gate $t_stop/$(length(phases))\n"; color = :yellow)
end

_ = camps_circuit_dynamics(ψ, gates, phases, χ_camps, obs, output; showprogress = true)

ev = cmps.expectation(ψ, obs)
append_datapoint(output, 0, real(ev))
_ = pauliprop_circuit_dynamics(ψ, 0, gates, phases, thl_pp, Nmax_pp, obs, output; showprogress = true)
return