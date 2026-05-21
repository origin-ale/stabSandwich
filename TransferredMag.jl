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

seed!(2)

N = 22
t = 11
ϕ = π/4
θ = π/4
μ = 100 # μ = 100 gives hard domain wall already, equivalent to μ = ∞
Nsamples = 50

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

layer_ends = layerends(N, t, xxz_circuit)

gates, phases = xxz_circuit(ϕ, θ, t, N)
gates, phases, layer_ends = dopeMagic(
  N, gates, phases, layer_ends, magic_syms, magic_inds, magic_prob)
  
printstyled("Running μ = $μ XXZ circuit dynamics until t = $t for \
N=$N, thl_pp = $thl_pp, Nmax_pp = $Nmax_pp.\n"; color = :cyan)

evs_pp = []
evs_cpp = []
prog = Progress(Nsamples; desc = "Computing ensemble averages…")
for it in 1:Nsamples
  ψ, onebitinds = domainwallstate(N, μ)
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

save_three_columns(layers, ["\n\n"], ["\n\n"], output)
save_three_columns(layers, ev_means_pp, ev_errs_pp, output)

return