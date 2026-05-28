using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!, MersenneTwister
using ProgressMeter
using Statistics

seed!(1)

N = 10
t = N ÷ 2
ϕ = π/4
θ = π/4
μ = 1 # μ = 100 gives hard domain wall already, equivalent to μ = ∞
Nsamples = 10

χ_campspp = 64
thl_campspp = 1e-15
Nmax_campspp = 100
χ_camps = χ_campspp
thl_pp = thl_campspp
Nmax_pp = 1_000_000
magic_prob = 1
magic_syms = [:Z for i = 1:N]
magic_inds = collect(1:N)
output = "output/TMDynamics.txt"

layer_ends = layerends(N, t, xxz_circuit)

gates, phases = xxz_circuit(ϕ, θ, t, N)
phases = subMagic(phases, magic_prob)

printstyled("Running μ = $μ XXZ circuit dynamics until t = $t for \
N=$N, thl_pp = $thl_pp, Nmax_pp = $Nmax_pp.\n"; color = :cyan)

evs_cpp = []
evs_pp = []
prog = Progress(Nsamples; desc = "Computing ensemble averages…")
for it in 1:Nsamples
  rng = MersenneTwister(it)
  ψ, onebitinds = domainwallstate(rng, N, μ)
  obs = transferredmagnetization(N, onebitinds)

  evs_it_cpp = campspp_circuit_dynamics(
    ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output;
    layer_ends = layer_ends)

  _, evs_it_pp = pauliprop_circuit_dynamics(
    onebitinds, gates, phases, thl_pp, Nmax_pp, obs, output; 
    layer_ends = layer_ends)

    global evs_cpp = isempty(evs_cpp) ? evs_it_cpp : hcat(evs_cpp, evs_it_cpp)
    global evs_pp = isempty(evs_pp) ? evs_it_pp : hcat(evs_pp, evs_it_pp)
    next!(prog; showvalues = [("ones on", onebitinds)])
end

ev_means_cpp = mean(evs_cpp, dims=2)
ev_errs_cpp = std(evs_cpp, dims=2)/sqrt(Nsamples)
ev_means_pp = mean(evs_pp, dims=2)
ev_errs_pp = std(evs_pp, dims=2)/sqrt(Nsamples)
layers = collect(0:length(layer_ends))

initialize_output(
  output, 
  "Transferred magnetization (avg over $Nsamples samples)", 
  Dict(
    "N" => N, 
    "μ" => μ,
    "Δ" => ϕ/θ,
    "p_dope" => magic_prob,
    "thl" => thl_pp,
    "Nmax" => Nmax_pp))

save_three_columns(layers, ev_means_cpp, ev_errs_cpp, output)
save_three_columns(["\n\n"], [""], [""], output)
save_three_columns(layers, ev_means_pp, ev_errs_pp, output)

return