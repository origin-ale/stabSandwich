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
using SHA
using Strided
using LinearAlgebra

Strided.disable_threads()
nthr=Threads.nthreads()

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

N = 14
t = N ÷ 2
ϕ = π/4
θ = π/4
μs = [0.3, 0.6, 1., 10.]
Nsamples = 30

χ_campspp = 128
thl_campspp = 1e-10
Nmax_campspp = 1000
magic_prob = 1
output = "output/TMD_$(N)_$(round(magic_prob))_$(Nsamples).txt"
output_full = "output/TMD_$(N)_$(round(magic_prob))_$(Nsamples)_full.txt"

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
N=$N, thl = $thl_campspp, Nmax = $Nmax_campspp.\nNsamples = $Nsamples, $nthr threads.\n"; color = :cyan)

prog = Progress(length(μs) * Nsamples; desc = "Computing…")

evs_by_μ = Vector{Any}(undef, length(μs))
sample_evs_by_μ = Vector{Any}(undef, length(μs))
sample_onebitinds_by_μ = Vector{Any}(undef, length(μs))

for μ_idx in eachindex(μs)
  μ = μs[μ_idx]
  sample_evs = Vector{Any}(undef, Nsamples)
  sample_onebitinds = Vector{Vector{Int}}(undef, Nsamples)

  @sync for it in 1:Nsamples
    Threads.@spawn begin
      rng = MersenneTwister(100_000 * μ_idx + it)
      ψ, onebitinds = domainwallstate(rng, N, μ)
      obs = transferredmagnetization(N, onebitinds)

      evs_it = campspp_circuit_dynamics(
        ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs;
        layer_ends = layer_ends)

      sample_evs[it] = evs_it
      sample_onebitinds[it] = copy(onebitinds)
      next!(prog)
    end
  end

  evs = hcat(sample_evs...)
  evs_by_μ[μ_idx] = evs
  sample_evs_by_μ[μ_idx] = sample_evs
  sample_onebitinds_by_μ[μ_idx] = sample_onebitinds
end

initialize_output(
  output_full,
  "Transferred magnetization (all $Nsamples samples)",
  Dict(
    "N" => N,
    "μs" => μs,
    "Δ" => ϕ/θ,
    "p_dope" => magic_prob,
    "χ" => χ_campspp,
    "thl" => thl_campspp,
    "Nmax" => Nmax_campspp))

open(output_full, "a") do full_io
  for (μ_idx, μ) in pairs(μs)
    println(full_io, "# μ = $μ")
    for (sample_idx, evs_it) in enumerate(sample_evs_by_μ[μ_idx])
      println(full_io, "# sample $sample_idx")
      for (layer_idx, ev) in enumerate(evs_it)
        println(full_io, "$(layer_idx - 1)\t$(ev)")
      end
      println(full_io)
    end
  end
end

for (μ_idx, μ) in pairs(μs)
  evs = evs_by_μ[μ_idx]
  ev_means = mean(evs, dims=2)
  ev_errs = std(evs, dims=2)/sqrt(Nsamples)
  layers = collect(0:length(layer_ends))

  save_three_columns(layers, ev_means, ev_errs, output)
  save_three_columns(["\n\n"], [""], [""], output)
end

return
