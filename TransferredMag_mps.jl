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


N = 16
t = N ÷ 2
ϕ = π/4
θ = π/4
μs = [0.3, 0.6, 1., 10.]
magic_probs = [1]
Nsamples = 10

dope_phase = 3π/16
dope_method = "on XY"

param_pairs = vec([(magic_prob, μ) for magic_prob in magic_probs, μ in μs])

thl = 1e-7

prefix = "mps_TMDxy"
output = "output/$(prefix)_$(N)_$(Nsamples).txt"
output_log = "output/$(prefix)_$(N)_$(Nsamples)_log.txt"
output_full = "output/$(prefix)_$(N)_$(Nsamples)_full.txt"
param_info = Dict(
  "N" => N,
  "Δ" => θ/ϕ,
  "doping" => "$(round(dope_phase; digits = 3)) $dope_method",
  "thl" => thl,
  "n.gates" => length(xxz_circuit(1,1, t, N)[1]))
obsname = "Transferred magnetization"
initialize_output(output, "$obsname (avg over $Nsamples samples)", param_info)
initialize_output(output_log, "$obsname", param_info)
initialize_output(output_full, "$obsname", param_info)

printstyled("Running MPS XXZ circuit dynamics over $(length(param_pairs)) (magic_prob, μ) pairs \
until t = $t for N=$N. thl = $thl\nNsamples = $Nsamples.\n"; color = :cyan)
prog = Progress(length(param_pairs) * Nsamples; desc = "Sampling…", enabled = true)
ProgressMeter.update!(prog, 0)


for pair_idx in eachindex(param_pairs)
  magic_prob, μ = param_pairs[pair_idx]
  sample_evs = Vector{Any}(undef, Nsamples)

  for it in 1:Nsamples
    rng = MersenneTwister(100_000 * pair_idx + it)

    layer_ends = layerends(N, t, xxz_circuit)
    gates, phases = xxz_circuit(ϕ, θ, t, N)
    phases = xy_magic(rng, phases, magic_prob; magicphase=dope_phase)

    ψ, onebitinds = domainwallstate(rng, N, μ)
    onebitinds_strs = fill("0", N)
    for i in onebitinds
      onebitinds_strs[i] = "1"
    end
    sites = siteinds("Qubit", N)
    ψ = MPS(sites, onebitinds_strs)
    obs = transferredmagnetization(N, onebitinds)

    _, evs_it = mps_circuit_dynamics(ψ, gates, phases, thl, obs, output_log; layer_ends=layer_ends)
    sample_evs[it] = evs_it
    next!(prog)
  end

  evs = stack_samples(sample_evs)
  save_full_samples(output_full, μ, magic_prob, evs)
  save_stats(output, evs, μ, magic_prob)
end

return
