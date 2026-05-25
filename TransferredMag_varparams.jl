using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!
using ProgressMeter
using Statistics

seed!(1)

N = 10
t = N ÷ 2
ϕ = π/4
θ = π/4
μs = [0.3, 0.6, 1, 10]
Nsamples = 25

χ_campspp = 64
thl_campspp = 1e-15
Nmax_campspp = 100
magic_prob = 1
output = "output/TMDynamics_varparams.txt"
output_full = "output/TMDynamics_varparams_full.txt"

layer_ends = layerends(N, t, xxz_circuit)

gates, phases = xxz_circuit(ϕ, θ, t, N)
phases = subMagic(phases, magic_prob)

initialize_output(
  output, 
  "Transferred magnetization (avg over $Nsamples samples)", 
  Dict(
    "N" => N, 
    "μs" => μs,
    "Δ" => ϕ/θ,
    "p_dope" => magic_prob,
    "χ" => χ_campspp,
    "thl" => thl_campspp,
    "Nmax" => Nmax_campspp))

printstyled("Running XXZ circuit dynamics until t = $t for \
N=$N, thl = $thl_campspp, Nmax_pp = $Nmax_campspp.\n"; color = :cyan)

prog = Progress(Nsamples*length(μs); desc = "Computing…")
for μ in μs
  evs = []
  for it in 1:Nsamples
    ψ, onebitinds = domainwallstate(N, μ)
    obs = transferredmagnetization(N, onebitinds)

    evs_it = campspp_circuit_dynamics(
      ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output_full;
      layer_ends = layer_ends)

    evs = isempty(evs) ? evs_it : hcat(evs, evs_it)
    next!(prog; showvalues = [("μ", μ),("sample", it)])
  end
  ev_means = mean(evs, dims=2)
  ev_errs = std(evs, dims=2)/sqrt(Nsamples)
  layers = collect(0:length(layer_ends))

  save_three_columns(layers, ev_means, ev_errs, output)
  save_three_columns(["\n\n"], [""], [""], output)
end

return