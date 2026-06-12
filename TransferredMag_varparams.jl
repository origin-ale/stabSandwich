using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Printf
using Random: seed!, MersenneTwister
using ProgressMeter
using SHA
using Strided
using LinearAlgebra

Strided.disable_threads()
nthr=Threads.nthreads()

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

N = 12
t = N ÷ 2
ϕ = π/4
θ = π/4
μs = [0.3, 0.6, 1., 10.]
Nsamples = 25
warn_on_prestop = true

χ_campspp = 128
thl_campspp = 1e-10
Nmax_campspp = 1000
magic_prob = 0
output = "output/TMD_$(N)_$(round(magic_prob))_$(Nsamples).txt"
output_log = "output/TMD_$(N)_$(round(magic_prob))_$(Nsamples)_log.txt"
output_full = "output/TMD_$(N)_$(round(magic_prob))_$(Nsamples)_full.txt"

layer_ends = layerends(N, t, xxz_circuit)

gates, phases = xxz_circuit(ϕ, θ, t, N)
dope_syms = fill(:Z, N)
dope_inds = collect(1:N)
dope_phase = π/8
# gates, phases, layer_ends = CampsPP.dopeMagic(N, gates, phases, layer_ends, dope_syms, dope_inds, magic_prob)
phases = subMagic(phases, magic_prob; magicphase = dope_phase)

param_info = Dict(
  "N" => N, 
  "μs" => μs,
  "Δ" => ϕ/θ,
  "p_dope" => magic_prob,
  "χ" => χ_campspp,
  "thl" => thl_campspp,
  "Nmax" => Nmax_campspp,
  "n.gates" => length(phases))
obsname = "Transferred magnetization"
initialize_output(output, "$obsname (avg over $Nsamples samples)", param_info)
initialize_output(output_log, "$obsname", param_info)
initialize_output(output_full, "$obsname", param_info)

printstyled("Running XXZ circuit dynamics with $magic_prob/1 doping until t = $t for \
N=$N, thl = $thl_campspp, Nmax = $Nmax_campspp.\nNsamples = $Nsamples, $nthr threads.\n"; color = :cyan)

prog = Progress(length(μs) * Nsamples; desc = "Computing…", enabled = true)

for μ_idx in eachindex(μs)
  μ = μs[μ_idx]
  sample_evs = Vector{Any}(undef, Nsamples)

  for it in 1:Nsamples
    rng = MersenneTwister(100_000 * μ_idx + it)
    ψ, onebitinds = domainwallstate(rng, N, μ)
    obs = transferredmagnetization(N, onebitinds)

    evs_it, t_stop = campspp_circuit_dynamics(
      ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, output_log;
      layer_ends = layer_ends)
    if warn_on_prestop && (t_stop < length(phases))
      printstyled("\rWARNING: sample $it for μ=$μ stopped at gate $t_stop/$(length(phases))    \n"; color = :yellow)
    end
    sample_evs[it] = evs_it
    next!(prog)
  end

  evs = stack_samples(sample_evs)
  save_full_samples(output_full, μ, evs)
  save_stats(output, evs)
end

return
