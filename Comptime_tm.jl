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
BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

# == Parameters ===============================================================
nmethods = 3 # CAMPS-PP, CAMPS, and MPS

magic_prob = 0
magic_mode = :xy # Dope on XX-YY with 3π/16 or on ZZ with π/3
μ = 0.6

Nmin = 2
spacing = 2
Nmax = 10

samples = 10

ϕ = π/4
θ = π/4

# CAMPS-PP (sandwich, ie. s-)
campspp_χ = 64 # s-CAMPS bond dimension that triggers s-PP
campspp_thl = 1e-10 # Threshold for truncation during s-PP
campspp_Pmax = 1000 # Maximum number of Paulis for s-PP

# CAMPS
camps_thl = 1e-7 # CAMPS SVD truncation threshold
camps_crit = :chi3 # CAMPS entanglement criterion (:chi3 or :entangle)
camps_strat = :snake # CAMPS disentangler strategy (:full, :brickwork, :snake)

# MPS
mps_thl = 1e-7 # MPS SVD truncation threshold

# output
output = "output/comptime_tm.txt"

# =============================================================================

magic_prob_str = @sprintf("%.2f", magic_prob)
Ns = Nmin:spacing:Nmax

param_info = Dict(
  "magic_prob" => magic_prob,
  "magic_mode" => magic_mode,
  "μ" => μ,
  "Nmin" => Nmin,
  "spacing" => spacing,
  "Nmax" => Nmax,
  "samples" => samples,
  "ϕ" => ϕ,
  "θ" => θ,
  "campspp_χ" => campspp_χ,
  "campspp_thl" => campspp_thl,
  "campspp_Pmax" => campspp_Pmax,
  "camps_thl" => camps_thl,
  "camps_crit" => camps_crit,
  "camps_strat" => camps_strat,
  "mps_thl" => mps_thl)

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

prefix = "output/comptime_p$(magic_prob_str)_"
suffix = "_log.txt"
campspp_log = prefix * "campspp" * suffix
campssrc_log = prefix * "campssrc" * suffix
mps_log = prefix * "mps" * suffix

times_methods = [Real[] for i in 1:nmethods]

function init_camps(N, seed)
  rng = MersenneTwister(seed)
  ψ, onebitinds = domainwallstate(rng, N, μ)
  tm = transferredmagnetization(N, onebitinds)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  phases = magic_doping(phases, magic_prob; magicphase = magic_phase)
  return ψ, tm, gates, phases
end

function init_mps(N, seed)
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
  phases = magic_doping(phases, magic_prob; magicphase = magic_phase)
  return ψ, tm, gates, phases
end

# == Main loop ================================================================
for N in Ns
  layer_ends = layerends(N, N/2, xxz_circuit)
  # == CAMPS-PP ===============================================================
  prog = Progress(samples; desc = "N=$N CAMPS-PP")

  ψ_wu, tm_wu, gates, phases = init_camps(N, 0)
  times_curr = Real[]
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
    ψ, tm, gates, phases = init_camps(N, i)
    (evs, tstop), time, _ = @timed campspp_circuit_dynamics(
      ψ, 
      campspp_χ, 
      campspp_thl, 
      campspp_Pmax, 
      gates, 
      phases, 
      tm,
      campspp_log,
      layer_ends = layer_ends)
    push!(evs_campspp, evs)
    push!(times_curr, time)
    next!(prog)
  end
  time_campspp = mean(times_curr)
  push!(times_methods[1], time_campspp)

  # == CAMPS ==================================================================
  prog = Progress(samples; desc = "N=$N CAMPS")

  ψ_wu, tm_wu, gates, phases = init_camps(N, 0)
  times_curr = Real[]
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
    ψ, tm, gates, phases = init_camps(N, i)
    (_, _, evs), time, _ = @timed campssrc_circuit_dynamics(
      ψ,
      gates,
      phases,
      camps_thl,
      tm,
      campssrc_log;
      criterion = camps_crit,
      strategy = camps_strat,
      layer_ends = layer_ends)
    push!(evs_camps, evs)
    push!(times_curr, time)
    next!(prog)
  end
  time_camps = mean(times_curr)
  push!(times_methods[2], time_camps)

  # == MPS ==================================================================
  prog = Progress(samples; desc = "N=$N MPS")

  ψ_wu, tm_wu, gates, phases = init_mps(N, 0)
  layer_ends = layerends(N, N ÷ 2, xxz_circuit)
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
    ψ, tm, gates, phases = init_mps(N, i)
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
  push!(times_methods[3], time_mps)

  if !all(isapprox(evs_campspp, evs_camps; atol = 1e-10))
    printstyled("WARNING: N = $N CAMPS-PP and CAMPS do not match\n";
      color=:yellow)
    println("-"^32)
    for idx in eachindex(evs_campspp)
      println("$(evs_campspp[idx])\n$(evs_camps[idx])\n")
    end
    println("-"^32)
  end
  if !all(isapprox(evs_campspp, evs_mps; atol = 1e-10))
    printstyled("WARNING: N = $N CAMPS-PP and MPS do not match\n";
      color=:yellow)
    println("-"^32)
    for idx in eachindex(evs_campspp)
      println("$(evs_campspp[idx])\n $(evs_mps[idx])\n")
    end
    println("-"^32)
  end

  initialize_output(output, "[ignore this row]", param_info)
  save_columns(output, Ns, times_methods...)
end

