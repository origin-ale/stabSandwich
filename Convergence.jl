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

N = 46
t = N ÷ 2
ϕ = π/4
θ = π/4
μ = 10

magic_prob = 0.007
dope_mode = "on XX-YY"
dope_phase = 3π/16

thl = 1
thl_step = 10
convergence = 1e-9
nsamples = 50

χ = 64
Nmax_pauli = 1_000
warn_on_prestop = true

output_log = "output/convergence_log.txt"
param_info = Dict(
  "N" => N,
  "Δ" => θ/ϕ,
  "doping" => "$(round(dope_phase; digits = 3)) $dope_mode",
  "χ" => χ,
  "Nmax_pauli" => Nmax_pauli,
  "n.gates" => length(xxz_circuit(1,1, t, N)[1]))
initialize_output(output_log, "Transferred magnetization", param_info)

printstyled(
  "Testing convergence of CAMPS-PP for XXZ Trotterization \
  with doping $(round(dope_phase; digits = 3)) $dope_mode \
  with p = $magic_prob until t = $t for N = $N.\n\
  χ = $χ, Nsamples = $nsamples, Nmax_pauli = $Nmax_pauli.\n";
  color = :cyan)

rng = MersenneTwister(0)

layer_ends = layerends(N, t, xxz_circuit)

ψ, onebitinds = domainwallstate(rng, N, μ)
obs = transferredmagnetization(N, onebitinds)

gates, phases = xxz_circuit(ϕ, θ, t, N)
phases = xy_magic(rng, phases, magic_prob; magicphase=dope_phase)

evs, _ = campspp_circuit_dynamics(
    ψ, χ, thl, Nmax_pauli, gates, phases, obs;
    layer_ends = layer_ends,
    showprogress=false)

last_evs = copy(evs)
fill!(last_evs, Inf)

diff = @. abs(last_evs - evs)

while maximum(diff) > convergence
  global thl /= thl_step
  global last_evs = copy(evs)

  open(output_log, "a") do f
    println(f, "========================= thl = $thl =========================\n")
  end

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

    global sample_evs[it], t_stop = campspp_circuit_dynamics(
      ψ, χ, thl, Nmax_pauli, gates, phases, obs, output_log,
      layer_ends = layer_ends)
    if warn_on_prestop && (t_stop < length(phases))
      printstyled("\rWARNING: sample $it stopped at gate \
      $t_stop/$(length(phases))\n"; color = :yellow)
    end
    next!(prog)
  end
  evs_its = stack_samples(sample_evs)
  global evs = [mean(skipmissing(row)) for row in eachrow(evs_its)]
  global last_evs = last_evs[1:length(evs)]
  global diff = @. abs(last_evs - evs)
  @show evs
  @show last_evs
  @printf "PP threshold %.1e -> max variation %.3e\n" thl maximum(diff)
end

outstr = @sprintf("PP converged at threshold %.1e, with max variation %.3e\n",
  thl, 
  maximum(diff))
printstyled(outstr, "\n"; color = :cyan)
return
