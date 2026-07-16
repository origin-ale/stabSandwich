using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Statistics
using Random: MersenneTwister
using ProgressMeter
using Printf

using Strided
using LinearAlgebra
Strided.disable_threads()
nthr = Threads.nthreads()
BLAS.set_num_threads(nthr)
ITensors.Strided.set_num_threads(nthr)

# == Parameters ===============================================================
# Methods to run, in output-column order. Comment out any you want to skip.
# Available: :campspp (CAMPS-PP), :camps (CAMPS), :mps (MPS), :pp (PP)
methods = [:campspp, :pp]

N = 46 # Fixed system size
magic_mode = :xy # Dope on XX-YY with 3π/16 or on ZZ with π/3
μ = 0.6

magic_prob_min = 0.000
magic_prob_spacing = 0.005
magic_prob_max = 0.050

samples = 25

ϕ = π/4
θ = π/4

# CAMPS-PP (sandwich, ie. s-)
campspp_χ = 128 # s-CAMPS bond dimension that triggers s-PP
campspp_thl = 1e-12 # Threshold for truncation during s-PP
campspp_Pmax = 10_000_000 # Maximum number of Paulis for s-PP

# CAMPS
camps_thl = 1e-12 # CAMPS SVD truncation threshold
camps_crit = :chi3 # CAMPS entanglement criterion (:chi3 or :entangle)
camps_strat = :snake # CAMPS disentangler strategy (:full, :brickwork, :snake)

# MPS
mps_thl = 1e-7 # MPS SVD truncation threshold

# Pauli propagation
pp_thl = 1e-10 # Threshold for truncation during PP
pp_Pmax = 10_000_000 # Maximum number of Paulis for PP

# output
output = "output/comptime_prob_tm.txt"
resources_prefix = "output/resources_prob_tm_"

# =============================================================================

all_methods = [:campspp, :camps, :mps, :pp]
let unknown = setdiff(methods, all_methods)
  isempty(unknown) || error("Unrecognised method(s): $unknown")
end
isempty(methods) && error("No methods selected to run")

magic_probs = magic_prob_min:magic_prob_spacing:magic_prob_max

param_info = Dict(
  "methods" => methods,
  "N" => N,
  "magic_mode" => magic_mode,
  "μ" => μ,
  "magic_prob_min" => magic_prob_min,
  "magic_prob_spacing" => magic_prob_spacing,
  "magic_prob_max" => magic_prob_max,
  "samples" => samples,
  "ϕ" => ϕ,
  "θ" => θ,
  "campspp_χ" => campspp_χ,
  "campspp_thl" => campspp_thl,
  "campspp_Pmax" => campspp_Pmax,
  "camps_thl" => camps_thl,
  "camps_crit" => camps_crit,
  "camps_strat" => camps_strat,
  "mps_thl" => mps_thl,
  "pp_thl" => pp_thl,
  "pp_Pmax" => pp_Pmax)

if magic_mode == :xy
  magic_phase = 3π/16
  magic_doping = xy_magic
  magic_txt = "3π/16 on XX-YY"
elseif magic_mode == :z
  magic_phase = 3π/16
  magic_doping = z_magic
  magic_txt = "π/3 on ZZ"
else
  error("Unrecognised magic_mode: $magic_mode")
end

prefix = "output/comptime_N$(N)_"
suffix = "_log.txt"
campspp_log = prefix * "campspp" * suffix
campssrc_log = prefix * "campssrc" * suffix
mps_log = prefix * "mps" * suffix
pp_log = prefix * "pp" * suffix

times_methods = Dict(m => Real[] for m in methods)
stds_methods = Dict(m => Real[] for m in methods)

# Per-layer bond dimension (:campspp, :camps) and n. of Paulis (:campspp, :pp),
# averaged over samples; one data block per probability
bd_outputs = Dict(m => resources_prefix * "bd_$m.txt"
  for m in methods if m in (:campspp, :camps))
np_outputs = Dict(m => resources_prefix * "NP_$m.txt"
  for m in methods if m in (:campspp, :pp))
for f in values(bd_outputs)
  initialize_output(f, "[layer, mean bond dim., std. err., max-sample bond dim.; one block per p]", param_info)
end
for f in values(np_outputs)
  initialize_output(f, "[layer, mean n. of Paulis, std. err., max-sample n. of Paulis; one block per p]", param_info)
end

function init_camps(N, seed, magic_prob)
  rng = MersenneTwister(seed)
  ψ, onebitinds = domainwallstate(rng, N, μ)
  tm = transferredmagnetization(N, onebitinds)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  doping_rng = MersenneTwister(hash((seed, :doping)))
  phases = magic_doping(doping_rng, phases, magic_prob; magicphase = magic_phase)
  return ψ, tm, gates, phases
end

function init_mps(N, seed, magic_prob)
  rng = MersenneTwister(seed)
  _, onebitinds = domainwallstate(rng, N, μ)
  tm = transferredmagnetization(N, onebitinds)

  onebitinds_strs = fill("0", N)
  for i in onebitinds
    onebitinds_strs[i] = "1"
  end
  sites = siteinds("Qubit", N)
  ψ = MPS(sites, onebitinds_strs)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  doping_rng = MersenneTwister(hash((seed, :doping)))
  phases = magic_doping(doping_rng, phases, magic_prob; magicphase = magic_phase)
  return ψ, tm, gates, phases
end

function init_pp(N, seed, magic_prob)
  rng = MersenneTwister(seed)
  _, onebitinds = domainwallstate(rng, N, μ)
  tm = transferredmagnetization(N, onebitinds)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  doping_rng = MersenneTwister(hash((seed, :doping)))
  phases = magic_doping(doping_rng, phases, magic_prob; magicphase = magic_phase)
  return onebitinds, tm, gates, phases
end

# == Main loop ================================================================
layer_ends = layerends(N, N ÷ 2, xxz_circuit)

for magic_prob in magic_probs
  magic_prob_str = @sprintf("%.3f", magic_prob)
  evs_campspp = evs_camps = evs_mps = evs_pp = nothing

  # == CAMPS-PP ===============================================================
  if :campspp in methods
    prog = Progress(samples; desc = "p=$magic_prob_str CAMPS-PP")

    ψ_wu, tm_wu, gates, phases = init_camps(N, 0, magic_prob)
    times_curr = Real[]
    bds_samples = Vector{Int}[]
    nps_samples = Vector{Int}[]
    evs_campspp = Vector{Real}[]
    _ = @timed campspp_circuit_dynamics(
      ψ_wu,
      campspp_χ,
      campspp_thl,
      campspp_Pmax,
      gates,
      phases,
      tm_wu,
      campspp_log;
      layer_ends = layer_ends)

    for i in 1:samples
      ψ, tm, gates, phases = init_camps(N, i, magic_prob)
      (evs, tstop, bds, nps), time, _ = @timed campspp_circuit_dynamics(
        ψ,
        campspp_χ,
        campspp_thl,
        campspp_Pmax,
        gates,
        phases,
        tm,
        campspp_log,
        layer_ends = layer_ends,
        track = true)
      if length(evs) < N/2+1
        printstyled("WARNING: CAMPS-PP sample $i for p=$magic_prob_str stopped early!")
      end
      push!(evs_campspp, evs)
      push!(times_curr, time)
      push!(bds_samples, bds)
      push!(nps_samples, nps)
      next!(prog)
    end
    time_campspp = mean(times_curr)
    std_campspp = std(times_curr) / sqrt(samples)
    push!(times_methods[:campspp], time_campspp)
    push!(stds_methods[:campspp], std_campspp)
    save_stats_maxcol(bd_outputs[:campspp], bds_samples, μ, magic_prob)
    save_stats_maxcol(np_outputs[:campspp], nps_samples, μ, magic_prob)
  end

  # == CAMPS ==================================================================
  if :camps in methods
    prog = Progress(samples; desc = "p=$magic_prob_str CAMPS")

    ψ_wu, tm_wu, gates, phases = init_camps(N, 0, magic_prob)
    times_curr = Real[]
    bds_samples = Vector{Int}[]
    evs_camps = Vector{Real}[]
    _ = @timed campssrc_circuit_dynamics(
      ψ_wu,
      gates,
      phases,
      camps_thl,
      tm_wu,
      campssrc_log;
      criterion = camps_crit,
      strategy = camps_strat,
      layer_ends = layer_ends)

    for i in 1:samples
      ψ, tm, gates, phases = init_camps(N, i, magic_prob)
      (_, _, evs, bds), time, _ = @timed campssrc_circuit_dynamics(
        ψ,
        gates,
        phases,
        camps_thl,
        tm,
        campssrc_log;
        criterion = camps_crit,
        strategy = camps_strat,
        layer_ends = layer_ends,
        track = true)
      push!(evs_camps, evs)
      push!(times_curr, time)
      push!(bds_samples, bds)
      next!(prog)
    end
    time_camps = mean(times_curr)
    std_camps = std(times_curr) / sqrt(samples)
    push!(times_methods[:camps], time_camps)
    push!(stds_methods[:camps], std_camps)
    save_stats_maxcol(bd_outputs[:camps], bds_samples, μ, magic_prob)
  end

  # == MPS ==================================================================
  if :mps in methods
    prog = Progress(samples; desc = "p=$magic_prob_str MPS")

    ψ_wu, tm_wu, gates, phases = init_mps(N, 0, magic_prob)
    times_curr = Real[]
    evs_mps = Vector{Real}[]
    _ = @timed mps_circuit_dynamics(
      ψ_wu,
      gates,
      phases,
      mps_thl,
      tm_wu,
      mps_log;
      layer_ends = layer_ends)

    for i in 1:samples
      ψ, tm, gates, phases = init_mps(N, i, magic_prob)
      (_, evs), time, _ = @timed mps_circuit_dynamics(
        ψ,
        gates,
        phases,
        mps_thl,
        tm,
        mps_log;
        layer_ends = layer_ends)
      push!(evs_mps, evs)
      push!(times_curr, time)
      next!(prog)
    end
    time_mps = mean(times_curr)
    std_mps = std(times_curr) / sqrt(samples)
    push!(times_methods[:mps], time_mps)
    push!(stds_methods[:mps], std_mps)
  end

  # == PP ====================================================================
  if :pp in methods
    prog = Progress(samples; desc = "p=$magic_prob_str PP")

    onebitinds_wu, tm_wu, gates, phases = init_pp(N, 0, magic_prob)
    times_curr = Real[]
    nps_samples = Vector{Int}[]
    evs_pp = Vector{Real}[]
    _ = @timed pauliprop_circuit_dynamics(
      onebitinds_wu,
      gates,
      phases,
      pp_thl,
      pp_Pmax,
      tm_wu,
      pp_log;
      layer_ends = layer_ends)

    for i in 1:samples
      onebitinds, tm, gates, phases = init_pp(N, i, magic_prob)
      (_, evs, nps), time, _ = @timed pauliprop_circuit_dynamics(
        onebitinds,
        gates,
        phases,
        pp_thl,
        pp_Pmax,
        tm,
        pp_log;
        layer_ends = layer_ends,
        track = true)
      if length(evs) < N/2+1
        printstyled("WARNING: PP sample $i for p=$magic_prob_str stopped early!")
      end
      push!(evs_pp, evs)
      push!(times_curr, time)
      push!(nps_samples, nps)
      next!(prog)
    end
    time_pp = mean(times_curr)
    std_pp = std(times_curr) / sqrt(samples)
    push!(times_methods[:pp], time_pp)
    push!(stds_methods[:pp], std_pp)
    save_stats_maxcol(np_outputs[:pp], nps_samples, μ, magic_prob)
  end

  # == Cross-checks (only between methods that were actually run) =============
  if !isnothing(evs_campspp) && !isnothing(evs_camps) &&
      !all(isapprox(evs_campspp, evs_camps; atol = 1e-6))
    printstyled("WARNING: p = $magic_prob_str CAMPS-PP and CAMPS do not match\n";
      color=:yellow)
    println("-"^32)
    for idx in eachindex(evs_campspp)
      println("$(evs_campspp[idx])\n$(evs_camps[idx])\n")
    end
    println("-"^32)
  end
  if !isnothing(evs_campspp) && !isnothing(evs_mps) &&
      !all(isapprox(evs_campspp, evs_mps; atol = 1e-6))
    printstyled("WARNING: p = $magic_prob_str CAMPS-PP and MPS do not match\n";
      color=:yellow)
    println("-"^32)
    for idx in eachindex(evs_campspp)
      println("$(evs_campspp[idx])\n $(evs_mps[idx])\n")
    end
    println("-"^32)
  end
  if !isnothing(evs_campspp) && !isnothing(evs_pp) &&
      !all(isapprox(evs_campspp, evs_pp; atol = 1e-6))
    printstyled("WARNING: p = $magic_prob_str CAMPS-PP and PP do not match\n";
      color=:yellow)
    println("-"^32)
    for idx in eachindex(evs_campspp)
      println("$(evs_campspp[idx])\n $(evs_pp[idx])\n")
    end
    println("-"^32)
  end

  initialize_output(output, "[ignore this row]", param_info)
  save_columns(output, magic_probs,
    (col for m in methods for col in (times_methods[m], stds_methods[m]))...)
end
