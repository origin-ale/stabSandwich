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


N = 14
t = N ÷ 2
ϕ = π/4
θ = π/4
Nsamples = 50

dope_phase = -π/8
phase_name = "pi8"

μ = 1
magic_prob = 0

χ = 128
thl = 1e-10
Nmax_pauli = 1000
warn_on_prestop = true

output = "output/TMD_expl_$(N)_$phase_name.txt"
param_info = Dict(
  "N" => N,
  "ϕ" => ϕ,
  "θ" => θ,
  "χ" => χ,
  "thl" => thl,
  "Nmax_pauli" => Nmax_pauli,
  "n.gates" => length(xxz_circuit(1,1, t, N)[1]),
  "doping phase" => round(dope_phase; digits=3))
obsname = "Transferred magnetization"
isfile(output) || initialize_output(output, "$obsname (avg over $Nsamples samples)", param_info)

printstyled("μ = $μ, p = $magic_prob, Nsamples = $Nsamples.\n"; color = :cyan)
prog = Progress(Nsamples; desc = "Sampling…", enabled = true)
ProgressMeter.update!(prog, 0)

sample_evs = Vector{Any}(undef, Nsamples)
for it in 1:Nsamples
  rng = MersenneTwister() # no explicit seed because I repeatedly run the same script

  layer_ends = layerends(N, t, xxz_circuit)
  gates, phases = xxz_circuit(ϕ, θ, t, N)
  phases = x_magic(rng, phases, magic_prob; magicphase=dope_phase)
  phases = y_magic(rng, phases, magic_prob; magicphase=dope_phase)

  ψ, onebitinds = domainwallstate(rng, N, μ)
  obs = transferredmagnetization(N, onebitinds)

  evs_it, t_stop = campspp_circuit_dynamics(
    ψ, χ, thl, Nmax_pauli, gates, phases, obs;
    layer_ends = layer_ends)
  if warn_on_prestop && (t_stop < length(phases))
    printstyled("\rWARNING: sample $it for μ=$μ stopped at gate $t_stop/$(length(phases))    \n"; color = :yellow)
  end
  sample_evs[it] = evs_it
  next!(prog)
end

evs = stack_samples(sample_evs)
append_stats(output, evs, μ, magic_prob)

return
