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

N = 12
t = N ÷ 2
ϕ = π/4
θ = π/4
μ = 10

magic_prob = .5
dope_mode = "on XX-YY"
dope_phase = 3π/16

thl = 1
thl_step = 10
convergence = 1e-9
nsamples = 50

Nmax_pauli = 1_000_000
warn_on_prestop = true

printstyled(
  "Testing convergence of PP for XXZ Trotterization \
  with doping $(round(dope_phase; digits = 3)) $dope_mode with p = $magic_prob \
  until t = $t for N = $N.\nNsamples = $nsamples, Nmax_pauli = $Nmax_pauli.\n";
  color = :cyan)

rng = MersenneTwister(0)

layer_ends = layerends(N, t, xxz_circuit)

ψ, onebitinds = domainwallstate(rng, N, μ)
obs = transferredmagnetization(N, onebitinds)

gates, phases = xxz_circuit(ϕ, θ, t, N)
phases = xy_magic(rng, phases, magic_prob; magicphase=dope_phase)

_, evs = pauliprop_circuit_dynamics(
    onebitinds,  gates, phases, thl, Nmax_pauli, obs;
    layer_ends = layer_ends,
    showprogress=false)

last_evs = copy(evs)
fill!(last_evs, Inf)

diff = @. abs(last_evs - evs)

while maximum(diff) > convergence
  global thl /= thl_step
  global last_evs = copy(evs)

  sample_evs = Vector{Any}(undef, nsamples)

  descstr = @sprintf "thl = %.1e…" thl
  prog = Progress(nsamples; desc = descstr, enabled = true)
  ProgressMeter.update!(prog, 0)

  for it in 1:nsamples
    local rng = MersenneTwister(it)

    local ψ, onebitinds = domainwallstate(rng, N, μ)
    local obs = transferredmagnetization(N, onebitinds)

    local gates, phases = xxz_circuit(ϕ, θ, t, N)
    local phases = xy_magic(rng, phases, magic_prob; magicphase=dope_phase)

    global t_stop, sample_evs[it] = pauliprop_circuit_dynamics(
      onebitinds,  gates, phases, thl, Nmax_pauli, obs,
      layer_ends = layer_ends)
    if warn_on_prestop && (t_stop < length(phases))
      printstyled("\rWARNING: sample $it stopped at gate \
      $t_stop/$(length(phases))\n"; color = :yellow)
    end
    next!(prog)
  end
  evs_its = stack_samples(sample_evs)
  global evs = [mean(skipmissing(row)) for row in eachrow(evs_its)]
  last_evs = last_evs[1:length(evs)]
  global diff = @. abs(last_evs - evs)
  @printf "PP threshold %.1e -> max variation %.3e\n" thl maximum(diff)
end

outstr = @sprintf("PP converged at threshold %.1e, with max variation %.3e\n",
  thl, 
  maximum(diff))
printstyled(outstr, "\n"; color = :cyan)
return
