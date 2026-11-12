using Revise
using CampsPP
using CampsPP.DisentangleCAMPS
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
# `methods`, `magic_probs` and `samples` can be overridden from the command
# line (see "Command-line overrides" below), e.g.
#   julia --project=. Comptime_prob_tm.jl methods=campspp p=0,0.01 samples=25

# Methods to run, in output-column order. Comment out any you want to skip.
# Available: :campspp (CAMPS-PP), :camps (CAMPS), :mps (MPS), :pp (PP)
methods = [:campspp]

N = 46 # Fixed system size
magic_mode = :xy # Dope on XX-YY with 3π/16 or on ZZ with π/3
μ = 0.6

# Doping probabilities, as an explicit list
magic_probs = [0.0, 0.01, 0.02, 0.04]

samples = 25
# Base RNG seed: sample i uses seed rng_seed + i (0 reproduces historical runs;
# change to draw an independent set of samples)
rng_seed = 100

ϕ = π/4
θ = π/4

# CAMPS-PP (sandwich, ie. s-)
campspp_χ = 64 # s-CAMPS bond dimension that triggers s-PP
campspp_thl = 1e-10 # Threshold for truncation during s-PP
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
gctime_output = "output/gctime_prob_tm.txt"
resources_prefix = "output/resources_prob_tm_"
# Mean expectation values; exponent of campspp_thl appended to file names
campspp_thl_exp = abs(round(Int, log10(campspp_thl)))
evs_prefix = "output/evs_prob_tm_"
# Doped rotation gates (Pauli string and qubits) per probability and sample
doped_output = "output/doped_prob_tm.txt"

# Per-sample initial states (one bitstring per sample, identical for every p
# and method) and CAMPS-PP switch points (one block per p)
initstates_output = "output/initstates_prob_tm.txt"
switch_output = "output/campspp_switch_prob_tm.txt"

# =============================================================================

# == Command-line overrides ===================================================
# Arguments are `key=value`; lists are comma-separated, without spaces:
#   methods=campspp,camps   p=0,0.01,0.02   samples=10
for arg in ARGS
  key, val = let kv = split(arg, "="; limit = 2)
    length(kv) == 2 || error("Argument `$arg` is not of the form key=value")
    kv
  end
  if key == "methods"
    global methods = Symbol.(split(val, ","))
  elseif key == "p"
    global magic_probs = parse.(Float64, split(val, ","))
  elseif key == "samples"
    global samples = parse(Int, val)
  else
    error("Unrecognised argument `$key` (expected methods, p or samples)")
  end
end

all_methods = [:campspp, :camps, :mps, :pp]
let unknown = setdiff(methods, all_methods)
  isempty(unknown) || error("Unrecognised method(s): $unknown")
end
isempty(methods) && error("No methods selected to run")
allunique(methods) || error("Repeated method(s) in $methods")
isempty(magic_probs) && error("No doping probabilities given")
all(p -> 0 <= p <= 1, magic_probs) || error("Probabilities must lie in [0, 1]: $magic_probs")
samples >= 1 || error("Need at least 1 sample, got $samples")

param_info = Dict(
  "methods" => join(methods, ","),
  "N" => N,
  "magic_mode" => magic_mode,
  "μ" => μ,
  "magic_probs" => join(magic_probs, ","),
  "samples" => samples,
  "rng_seed" => rng_seed,
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
  magic_phase = π/3
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
gctimes_methods = Dict(m => Real[] for m in methods)
gcstds_methods = Dict(m => Real[] for m in methods)

# Quartiles and extrema over samples, alongside the means above;
# q1/median/q3/min/max ordered for gnuplot `with candlesticks`
timeq1_methods = Dict(m => Real[] for m in methods)
timemed_methods = Dict(m => Real[] for m in methods)
timeq3_methods = Dict(m => Real[] for m in methods)
timemin_methods = Dict(m => Real[] for m in methods)
timemax_methods = Dict(m => Real[] for m in methods)
gcq1_methods = Dict(m => Real[] for m in methods)
gcmed_methods = Dict(m => Real[] for m in methods)
gcq3_methods = Dict(m => Real[] for m in methods)
gcmin_methods = Dict(m => Real[] for m in methods)
gcmax_methods = Dict(m => Real[] for m in methods)

""" Push the quartiles and extrema of `vals` onto the matching per-method
columns. """
function push_quartiles!(q1s, meds, q3s, mins, maxs, m, vals)
  q1, med, q3 = quantiles(vals)
  push!(q1s[m], q1)
  push!(meds[m], med)
  push!(q3s[m], q3)
  push!(mins[m], minimum(vals))
  push!(maxs[m], maximum(vals))
end

""" Indices of the samples whose ⟨O⟩ trajectories in `as` and `bs` differ by \
more than `atol` at some layer, or have different lengths (one method stopped \
early). """
mismatched_samples(as, bs; atol) = [i for i in eachindex(as, bs)
  if length(as[i]) != length(bs[i]) || !all(isapprox.(as[i], bs[i]; atol = atol))]

# Per-layer bond dimension (:campspp, :camps) and n. of Paulis (:campspp, :pp),
# averaged over samples; one data block per probability
bd_outputs = Dict(m => resources_prefix * "bd_$m.txt"
  for m in methods if m in (:campspp, :camps))
np_outputs = Dict(m => resources_prefix * "NP_$m.txt"
  for m in methods if m in (:campspp, :pp))
pw_outputs = Dict(m => resources_prefix * "PW_$m.txt"
  for m in methods if m in (:campspp, :pp))
pc_outputs = Dict(m => resources_prefix * "PC_$m.txt"
  for m in methods if m in (:campspp, :pp))
for f in values(bd_outputs)
  initialize_output(f, "[layer, mean bond dim., std. err., max-sample bond dim., n. of samples, q1, median, q3, min, max; one block per p]", param_info)
end
for f in values(np_outputs)
  initialize_output(f, "[layer, mean n. of Paulis, std. err., max-sample n. of Paulis, n. of samples, q1, median, q3, min, max; one block per p]", param_info)
end
for f in values(pw_outputs)
  initialize_output(f, "[layer, mean avg. Pauli weight, std. err., max-sample avg. Pauli weight, n. of samples, q1, median, q3, min, max; one block per p]", param_info)
end
for f in values(pc_outputs)
  initialize_output(f, "[layer, mean avg. |Pauli coeff.|, std. err., max-sample avg. |Pauli coeff.|, n. of samples, q1, median, q3, min, max; one block per p]", param_info)
end

# Unaveraged per-sample resources: layer-resolved, one block per sample
bd_full_outputs = Dict(m => resources_prefix * "bd_full_$m.txt"
  for m in keys(bd_outputs))
np_full_outputs = Dict(m => resources_prefix * "NP_full_$m.txt"
  for m in keys(np_outputs))
pw_full_outputs = Dict(m => resources_prefix * "PW_full_$m.txt"
  for m in keys(pw_outputs))
pc_full_outputs = Dict(m => resources_prefix * "PC_full_$m.txt"
  for m in keys(pc_outputs))
for (f, name) in Iterators.flatten((
    ((f, "bond dim.") for f in values(bd_full_outputs)),
    ((f, "n. of Paulis") for f in values(np_full_outputs)),
    ((f, "avg. Pauli weight") for f in values(pw_full_outputs)),
    ((f, "avg. |Pauli coeff.|") for f in values(pc_full_outputs))))
  initialize_output(f, "[layer, $name; one block per sample, one set per p]", param_info)
end

# Peak resource reached by each sample over its own evolution
bd_peak_outputs = Dict(m => resources_prefix * "bd_peak_$m.txt"
  for m in keys(bd_outputs))
np_peak_outputs = Dict(m => resources_prefix * "NP_peak_$m.txt"
  for m in keys(np_outputs))
pw_peak_outputs = Dict(m => resources_prefix * "PW_peak_$m.txt"
  for m in keys(pw_outputs))
pc_peak_outputs = Dict(m => resources_prefix * "PC_peak_$m.txt"
  for m in keys(pc_outputs))
for (f, name) in Iterators.flatten((
    ((f, "peak bond dim.") for f in values(bd_peak_outputs)),
    ((f, "peak n. of Paulis") for f in values(np_peak_outputs)),
    ((f, "peak avg. Pauli weight") for f in values(pw_peak_outputs)),
    ((f, "peak avg. |Pauli coeff.|") for f in values(pc_peak_outputs))))
  initialize_output(f, "[sample, $name; one block per p]", param_info)
end

# Unaveraged per-sample wall-clock and GC times
times_full_outputs = Dict(m => "output/comptime_full_prob_tm_$m.txt" for m in methods)
for f in values(times_full_outputs)
  initialize_output(f, "[sample, time (s), GC time (s); one block per p]", param_info)
end

# Per-layer mean expectation value, averaged over samples; one block per p
evs_outputs = Dict(m => evs_prefix * "$(m)_$(campspp_thl_exp).txt" for m in methods)
for f in values(evs_outputs)
  initialize_output(f, "[layer, mean expectation value, std. err.; one block per p]", param_info)
end

# Unaveraged per-sample expectation values, one block per sample
evs_full_outputs = Dict(m => evs_prefix * "$(m)_$(campspp_thl_exp)_full.txt"
  for m in methods)
for f in values(evs_full_outputs)
  initialize_output(f, "[layer, expectation value; one block per sample, one set per p]", param_info)
end

initialize_output(doped_output,
  "[doped gates: count line, then one gate position + cycle + Pauli string + qubits per gate; one set per sample, one block per p]",
  param_info)

if :campspp in methods
  n_gates = length(first(xxz_circuit(ϕ, θ, N/2, N)))
  initialize_output(switch_output,
    "[sample, switch gate (last gate applied in CAMPS; PP starts at the next one; $n_gates = never switched), n. of magic qubits k at the switch, PP stop gate ($n_gates = end of circuit); one block per p]",
    param_info)
end

# Initial states are drawn once per sample from seed `rng_seed + i` (index 0 is
# the warm-up sample) and saved, so a replay can rebuild them from the file
# instead of from the RNG. They do not depend on p or on the method.
sample_onebitinds(i) = domainwallstate(MersenneTwister(rng_seed + i), N, μ)[2]
onebitinds_wu = sample_onebitinds(0)
onebitinds_samples = [sample_onebitinds(i) for i in 1:samples]

initialize_output(initstates_output,
  "[sample, initial bitstring (character j is qubit j; 1 = excited), indices of excited qubits]",
  param_info)
open(initstates_output, "a") do io
  for (i, inds) in enumerate(onebitinds_samples)
    bits = join(j in inds ? '1' : '0' for j in 1:N)
    println(io, i, " ", bits, " ", join(inds, " "))
  end
end

# The doping `mask` is generated once per (sample, probability) in the main loop
# and passed in, so the doped gates are identical across methods for a sample.
function init_camps(N, onebitinds, mask)
  ψ = computationalcamps(N, onebitinds)
  tm = transferredmagnetization(N, onebitinds)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  phases = magic_doping(phases, mask; magicphase = magic_phase)
  return ψ, tm, gates, phases
end

function init_mps(N, onebitinds, mask)
  tm = transferredmagnetization(N, onebitinds)

  onebitinds_strs = fill("0", N)
  for i in onebitinds
    onebitinds_strs[i] = "1"
  end
  sites = siteinds("Qubit", N)
  ψ = MPS(sites, onebitinds_strs)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  phases = magic_doping(phases, mask; magicphase = magic_phase)
  return ψ, tm, gates, phases
end

function init_pp(N, onebitinds, mask)
  tm = transferredmagnetization(N, onebitinds)

  gates, phases = xxz_circuit(ϕ, θ, N/2, N)
  phases = magic_doping(phases, mask; magicphase = magic_phase)
  return onebitinds, tm, gates, phases
end

# == Main loop ================================================================
layer_ends = layerends(N, N ÷ 2, xxz_circuit)

for magic_prob in magic_probs
  magic_prob_str = @sprintf("%.3f", magic_prob)
  evs_campspp = evs_camps = evs_mps = evs_pp = nothing

  # Doping masks are fixed by the sample seed *and* probability, so they are
  # identical across methods for a given sample yet independent across
  # probabilities (no monotonic nesting). Index 0 is the warm-up sample.
  gates_ref, phases_ref = xxz_circuit(ϕ, θ, N/2, N)
  sample_mask(i) = doping_mask(MersenneTwister(hash((rng_seed + i, magic_prob))),
    phases_ref, magic_prob)
  masks_samples = [sample_mask(i) for i in 1:samples]
  mask_wu = sample_mask(0)

  # Record the resulting doped gates once per probability
  doped_inds_samples = map(masks_samples) do mask
    doped_phases = magic_doping(phases_ref, mask; magicphase = magic_phase)
    findall(doped_phases .!= phases_ref)
  end
  save_doped_gates(doped_output, gates_ref, doped_inds_samples, layer_ends, μ,
    magic_prob)

  # == CAMPS-PP ===============================================================
  if :campspp in methods
    prog = Progress(samples; desc = "p=$magic_prob_str CAMPS-PP")

    ψ_wu, tm_wu, gates, phases = init_camps(N, onebitinds_wu, mask_wu)
    times_curr = Real[]
    gctimes_curr = Real[]
    bds_samples = Vector{Int}[]
    nps_samples = Vector{Union{Missing, Int}}[]
    pws_samples = Vector{Union{Missing, Float64}}[]
    pcs_samples = Vector{Union{Missing, Float64}}[]
    evs_campspp = Vector{Real}[]
    switch_gates = Int[]
    switch_ks = Int[]
    stop_gates = Int[]
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
      ψ, tm, gates, phases = init_camps(N, onebitinds_samples[i], masks_samples[i])
      (evs, tstop, bds, nps, pws, pcs, s_switch, k_switch), time, _, gctime = @timed campspp_circuit_dynamics(
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
      push!(gctimes_curr, gctime)
      push!(bds_samples, bds)
      push!(switch_gates, s_switch)
      push!(switch_ks, k_switch)
      push!(stop_gates, tstop)
      # PP starts at the (sample-dependent) switch layer, one past the last
      # CAMPS-recorded layer; pad so that row k is circuit layer k
      push!(nps_samples, [fill(missing, length(bds) - 1); nps])
      push!(pws_samples, [fill(missing, length(bds) - 1); pws])
      push!(pcs_samples, [fill(missing, length(bds) - 1); pcs])
      next!(prog)
    end
    time_campspp = mean(times_curr)
    std_campspp = std(times_curr) / sqrt(samples)
    push!(times_methods[:campspp], time_campspp)
    push!(stds_methods[:campspp], std_campspp)
    push!(gctimes_methods[:campspp], mean(gctimes_curr))
    push!(gcstds_methods[:campspp], std(gctimes_curr) / sqrt(samples))
    push_quartiles!(timeq1_methods, timemed_methods, timeq3_methods,
      timemin_methods, timemax_methods, :campspp, times_curr)
    push_quartiles!(gcq1_methods, gcmed_methods, gcq3_methods,
      gcmin_methods, gcmax_methods, :campspp, gctimes_curr)
    save_stats_maxcol(bd_outputs[:campspp], bds_samples, μ, magic_prob)
    save_stats_maxcol(np_outputs[:campspp], nps_samples, μ, magic_prob)
    save_stats_maxcol(pw_outputs[:campspp], pws_samples, μ, magic_prob)
    save_stats_maxcol(pc_outputs[:campspp], pcs_samples, μ, magic_prob)
    save_stats(evs_outputs[:campspp], stack_samples(evs_campspp), μ, magic_prob)
    save_per_sample(times_full_outputs[:campspp], μ, magic_prob, times_curr, gctimes_curr)
    save_per_sample(switch_output, μ, magic_prob, switch_gates, switch_ks, stop_gates)
    save_full_samples(bd_full_outputs[:campspp], μ, magic_prob, stack_samples(bds_samples))
    save_full_samples(np_full_outputs[:campspp], μ, magic_prob, stack_samples(nps_samples))
    save_full_samples(pw_full_outputs[:campspp], μ, magic_prob, stack_samples(pws_samples))
    save_full_samples(pc_full_outputs[:campspp], μ, magic_prob, stack_samples(pcs_samples))
    save_full_samples(evs_full_outputs[:campspp], μ, magic_prob, stack_samples(evs_campspp))
    save_peaks(bd_peak_outputs[:campspp], bds_samples, μ, magic_prob)
    save_peaks(np_peak_outputs[:campspp], nps_samples, μ, magic_prob)
    save_peaks(pw_peak_outputs[:campspp], pws_samples, μ, magic_prob)
    save_peaks(pc_peak_outputs[:campspp], pcs_samples, μ, magic_prob)
  end

  # == CAMPS ==================================================================
  if :camps in methods
    prog = Progress(samples; desc = "p=$magic_prob_str CAMPS")

    ψ_wu, tm_wu, gates, phases = init_camps(N, onebitinds_wu, mask_wu)
    times_curr = Real[]
    gctimes_curr = Real[]
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
      ψ, tm, gates, phases = init_camps(N, onebitinds_samples[i], masks_samples[i])
      (_, _, evs, bds), time, _, gctime = @timed campssrc_circuit_dynamics(
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
      push!(gctimes_curr, gctime)
      push!(bds_samples, bds)
      next!(prog)
    end
    time_camps = mean(times_curr)
    std_camps = std(times_curr) / sqrt(samples)
    push!(times_methods[:camps], time_camps)
    push!(stds_methods[:camps], std_camps)
    push!(gctimes_methods[:camps], mean(gctimes_curr))
    push!(gcstds_methods[:camps], std(gctimes_curr) / sqrt(samples))
    push_quartiles!(timeq1_methods, timemed_methods, timeq3_methods,
      timemin_methods, timemax_methods, :camps, times_curr)
    push_quartiles!(gcq1_methods, gcmed_methods, gcq3_methods,
      gcmin_methods, gcmax_methods, :camps, gctimes_curr)
    save_stats_maxcol(bd_outputs[:camps], bds_samples, μ, magic_prob)
    save_stats(evs_outputs[:camps], stack_samples(evs_camps), μ, magic_prob)
    save_per_sample(times_full_outputs[:camps], μ, magic_prob, times_curr, gctimes_curr)
    save_full_samples(bd_full_outputs[:camps], μ, magic_prob, stack_samples(bds_samples))
    save_full_samples(evs_full_outputs[:camps], μ, magic_prob, stack_samples(evs_camps))
    save_peaks(bd_peak_outputs[:camps], bds_samples, μ, magic_prob)
  end

  # == MPS ==================================================================
  if :mps in methods
    prog = Progress(samples; desc = "p=$magic_prob_str MPS")

    ψ_wu, tm_wu, gates, phases = init_mps(N, onebitinds_wu, mask_wu)
    times_curr = Real[]
    gctimes_curr = Real[]
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
      ψ, tm, gates, phases = init_mps(N, onebitinds_samples[i], masks_samples[i])
      (_, evs), time, _, gctime = @timed mps_circuit_dynamics(
        ψ,
        gates,
        phases,
        mps_thl,
        tm,
        mps_log;
        layer_ends = layer_ends)
      push!(evs_mps, evs)
      push!(times_curr, time)
      push!(gctimes_curr, gctime)
      next!(prog)
    end
    time_mps = mean(times_curr)
    std_mps = std(times_curr) / sqrt(samples)
    push!(times_methods[:mps], time_mps)
    push!(stds_methods[:mps], std_mps)
    push!(gctimes_methods[:mps], mean(gctimes_curr))
    push!(gcstds_methods[:mps], std(gctimes_curr) / sqrt(samples))
    push_quartiles!(timeq1_methods, timemed_methods, timeq3_methods,
      timemin_methods, timemax_methods, :mps, times_curr)
    push_quartiles!(gcq1_methods, gcmed_methods, gcq3_methods,
      gcmin_methods, gcmax_methods, :mps, gctimes_curr)
    save_stats(evs_outputs[:mps], stack_samples(evs_mps), μ, magic_prob)
    save_per_sample(times_full_outputs[:mps], μ, magic_prob, times_curr, gctimes_curr)
    save_full_samples(evs_full_outputs[:mps], μ, magic_prob, stack_samples(evs_mps))
  end

  # == PP ====================================================================
  if :pp in methods
    prog = Progress(samples; desc = "p=$magic_prob_str PP")

    obi_wu, tm_wu, gates, phases = init_pp(N, onebitinds_wu, mask_wu)
    times_curr = Real[]
    gctimes_curr = Real[]
    nps_samples = Vector{Int}[]
    pws_samples = Vector{Float64}[]
    pcs_samples = Vector{Float64}[]
    evs_pp = Vector{Real}[]
    _ = @timed pauliprop_circuit_dynamics(
      obi_wu,
      gates,
      phases,
      pp_thl,
      pp_Pmax,
      tm_wu,
      pp_log;
      layer_ends = layer_ends)

    for i in 1:samples
      onebitinds, tm, gates, phases = init_pp(N, onebitinds_samples[i], masks_samples[i])
      (_, evs, nps, pws, pcs), time, _, gctime = @timed pauliprop_circuit_dynamics(
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
      push!(gctimes_curr, gctime)
      push!(nps_samples, nps)
      push!(pws_samples, pws)
      push!(pcs_samples, pcs)
      next!(prog)
    end
    time_pp = mean(times_curr)
    std_pp = std(times_curr) / sqrt(samples)
    push!(times_methods[:pp], time_pp)
    push!(stds_methods[:pp], std_pp)
    push!(gctimes_methods[:pp], mean(gctimes_curr))
    push!(gcstds_methods[:pp], std(gctimes_curr) / sqrt(samples))
    push_quartiles!(timeq1_methods, timemed_methods, timeq3_methods,
      timemin_methods, timemax_methods, :pp, times_curr)
    push_quartiles!(gcq1_methods, gcmed_methods, gcq3_methods,
      gcmin_methods, gcmax_methods, :pp, gctimes_curr)
    save_stats_maxcol(np_outputs[:pp], nps_samples, μ, magic_prob)
    save_stats_maxcol(pw_outputs[:pp], pws_samples, μ, magic_prob)
    save_stats_maxcol(pc_outputs[:pp], pcs_samples, μ, magic_prob)
    save_stats(evs_outputs[:pp], stack_samples(evs_pp), μ, magic_prob)
    save_per_sample(times_full_outputs[:pp], μ, magic_prob, times_curr, gctimes_curr)
    save_full_samples(np_full_outputs[:pp], μ, magic_prob, stack_samples(nps_samples))
    save_full_samples(pw_full_outputs[:pp], μ, magic_prob, stack_samples(pws_samples))
    save_full_samples(pc_full_outputs[:pp], μ, magic_prob, stack_samples(pcs_samples))
    save_full_samples(evs_full_outputs[:pp], μ, magic_prob, stack_samples(evs_pp))
    save_peaks(np_peak_outputs[:pp], nps_samples, μ, magic_prob)
    save_peaks(pw_peak_outputs[:pp], pws_samples, μ, magic_prob)
    save_peaks(pc_peak_outputs[:pp], pcs_samples, μ, magic_prob)
  end

  # == Cross-checks (only between methods that were actually run) =============
  # Per sample and per layer; a sample where one method stopped early (ragged
  # trajectories) counts as a mismatch instead of throwing DimensionMismatch
  for (name, evs_other) in (("CAMPS", evs_camps), ("MPS", evs_mps), ("PP", evs_pp))
    (isnothing(evs_campspp) || isnothing(evs_other)) && continue
    bad = mismatched_samples(evs_campspp, evs_other; atol = 1e-6)
    isempty(bad) && continue
    printstyled("WARNING: p = $magic_prob_str CAMPS-PP and $name do not match (samples $bad)\n";
      color=:yellow)
    println("-"^32)
    for idx in bad
      println("$(evs_campspp[idx])\n$(evs_other[idx])\n")
    end
    println("-"^32)
  end

  # Seven columns per method: mean, std. err., then q1/median/q3/min/max
  # (ordered for gnuplot `with candlesticks`, whiskers at the extrema)
  initialize_output(output,
    "[p, then per method: mean time (s), std. err., q1, median, q3, min, max]",
    param_info)
  save_columns(output, magic_probs,
    (col for m in methods for col in (times_methods[m], stds_methods[m],
      timeq1_methods[m], timemed_methods[m], timeq3_methods[m],
      timemin_methods[m], timemax_methods[m]))...)

  initialize_output(gctime_output,
    "[p, then per method: mean GC time (s), std. err., q1, median, q3, min, max]",
    param_info)
  save_columns(gctime_output, magic_probs,
    (col for m in methods for col in (gctimes_methods[m], gcstds_methods[m],
      gcq1_methods[m], gcmed_methods[m], gcq3_methods[m],
      gcmin_methods[m], gcmax_methods[m]))...)
end
