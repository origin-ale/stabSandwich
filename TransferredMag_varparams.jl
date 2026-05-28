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
using Strided
using LinearAlgebra

Strided.disable_threads()
nthr=Threads.nthreads()

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

seed!(1)

N = 22
t = N ÷ 2
ϕ = π/4
θ = π/4
μs = [0.3, 0.6, 1, 10]
Nsamples = 25

χ_campspp = 128
thl_campspp = 1e-10
Nmax_campspp = 1_000_000_000
magic_prob = 0
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

prog = Progress(length(μs)*Nsamples; desc = "Computing…")
printstyled("Running XXZ circuit dynamics until t = $t for \
N=$N, thl = $thl_campspp, Nmax = $Nmax_campspp.\nNsamples = $Nsamples, $nthr threads.\n"; color = :cyan)

evs_by_μ = Vector{Any}(undef, length(μs))
full_outputs = Vector{String}(undef, length(μs))

for μ_idx in eachindex(μs)
  μ = μs[μ_idx]
  temp_output_full = tempname()
  full_outputs[μ_idx] = temp_output_full
  sample_outputs = Vector{String}(undef, Nsamples)
  sample_evs = Vector{Any}(undef, Nsamples)

  Threads.@threads for it in 1:Nsamples
    ψ, onebitinds = domainwallstate(N, μ)
    obs = transferredmagnetization(N, onebitinds)

    temp_output = tempname()
    sample_outputs[it] = temp_output

    evs_it = campspp_circuit_dynamics(
      ψ, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, obs, temp_output;
      layer_ends = layer_ends)

    sample_evs[it] = evs_it
    next!(prog)
  end

  open(temp_output_full, "w") do full_io
    for temp_output in sample_outputs
      isfile(temp_output) || continue
      open(temp_output, "r") do temp_io
        write(full_io, read(temp_io))
      end
      rm(temp_output; force = true)
    end
  end

  evs = hcat(sample_evs...)
  evs_by_μ[μ_idx] = evs
end

open(output_full, "w") do full_io
  for temp_output in full_outputs
    isfile(temp_output) || continue
    open(temp_output, "r") do temp_io
      write(full_io, read(temp_io))
    end
    rm(temp_output; force = true)
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
