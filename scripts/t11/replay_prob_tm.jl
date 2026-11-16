# T11: replay the Phase 0 `Comptime_prob_tm.jl` samples on the new (MagicTensors
# `QubitCAMPS`) stack, sample by sample, for the old-vs-new paired comparison
# (`compare_prob_tm.jl`). Run from the repository root:
#
#   julia --project=. scripts/t11/replay_prob_tm.jl run=output/runs/<id> \
#     [methods=campspp] [p=0,0.01] [samples=1,2,3] [inject=none|phase|mask] [out=<dir>]
#
#   run      Phase 0 run directory (required). Its inputs are read from
#            <run>/files/output/ and checked against <run>/SHA256SUMS first.
#   methods  campspp and/or camps (default: the run's own `methods`)
#   p        doping probabilities (default: the run's own `magic_probs`)
#   samples  sample ids, 1-based as in the run (default: all of the run's samples)
#   inject   sensitivity check only: `phase` multiplies the magic phase by 1.05,
#            `mask` moves every doped gate pair one XXZ triple later. Default none.
#   out      output directory (default output/t11/<run id>[_inject-<inject>])
#
# Nothing is drawn from an RNG: each sample is rebuilt from its saved initial bitstring
# (`initstates_prob_tm.txt`) and its saved doped gates (`doped_prob_tm.txt`). Every
# parameter comes from the run's header line, not from this file. The outputs use the
# Phase 0 file names and the same `Utility.jl` writers, so one parser reads old and new;
# the header adds `sample_ids`, `replay_of`, `inject` and `stack`. The replay also writes
# the initial states and doped gates it actually used, so the comparison can check that
# the circuits are the same.

using CampsPP
import MagicTensors as mt
import PauliPropagation as pp
using ITensors, ITensorMPS

using Dates
using Printf
using Random: MersenneTwister
using SHA

using Strided
using LinearAlgebra
# Same threading set-up as Comptime_prob_tm.jl, so the times are comparable
Strided.disable_threads()
nthr = Threads.nthreads()
BLAS.set_num_threads(nthr)
ITensors.Strided.set_num_threads(nthr)

logline(msg) = (println("[", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "] ", msg); flush(stdout))

# == Arguments ================================================================
args = Dict{String, String}()
for arg in ARGS
  kv = split(arg, "="; limit = 2)
  length(kv) == 2 || error("Argument `$arg` is not of the form key=value")
  kv[1] in ("run", "methods", "p", "samples", "inject", "out") ||
    error("Unrecognised argument `$(kv[1])` (expected run, methods, p, samples, inject, out)")
  haskey(args, kv[1]) && error("Argument `$(kv[1])` given twice")
  args[kv[1]] = kv[2]
end
haskey(args, "run") || error("Missing run=<Phase 0 run directory>")
rundir = rstrip(args["run"], '/')
runid = basename(rundir)
filesdir = joinpath(rundir, "files", "output")
inject = Symbol(get(args, "inject", "none"))
inject in (:none, :phase, :mask) || error("inject must be none, phase or mask, got $inject")

# == Input integrity ==========================================================
initstates_file = joinpath(filesdir, "initstates_prob_tm.txt")
doped_file = joinpath(filesdir, "doped_prob_tm.txt")

let sums = Dict{String, String}()
  for line in eachline(joinpath(rundir, "SHA256SUMS"))
    h, f = split(line; limit = 2)
    sums[strip(f)] = h
  end
  for f in (initstates_file, doped_file)
    rel = relpath(f, rundir)
    haskey(sums, rel) || error("$rel is not listed in $rundir/SHA256SUMS")
    h = bytes2hex(open(sha256, f))
    h == sums[rel] || error("Checksum mismatch for $f")
  end
end

# == Parameters from the run's header =========================================
""" Parse the `key=value` header line written by `initialize_output`. """
function read_header(file)
  line = first(eachline(file))
  startswith(line, "#") || error("$file has no header line")
  return Dict(String(kv[1]) => String(kv[2])
    for kv in (split(t, "="; limit = 2) for t in split(line[2:end])) if length(kv) == 2)
end

hdr = read_header(initstates_file)
read_header(doped_file) == hdr || error("Header lines of the two input files differ")

N = parse(Int, hdr["N"])
μ = parse(Float64, hdr["μ"])
ϕ = parse(Float64, hdr["ϕ"])
θ = parse(Float64, hdr["θ"])
rng_seed = parse(Int, hdr["rng_seed"])
nsamples_run = parse(Int, hdr["samples"])
probs_run = parse.(Float64, split(hdr["magic_probs"], ","))
methods_run = Symbol.(split(hdr["methods"], ","))
magic_mode = Symbol(hdr["magic_mode"])
campspp_χ = parse(Int, hdr["campspp_χ"])
campspp_thl = parse(Float64, hdr["campspp_thl"])
campspp_Pmax = parse(Int, hdr["campspp_Pmax"])
camps_thl = parse(Float64, hdr["camps_thl"])
camps_crit = Symbol(hdr["camps_crit"])
camps_strat = Symbol(hdr["camps_strat"])

# The magic phase isn't in the header; Comptime_prob_tm.jl fixes it per mode
# (commit 2c4bf8e, lines 128-138). Only :xy has Phase 0 runs.
magic_mode == :xy || error("Only :xy doping is replayed (Phase 0 has no :z runs)")
magic_phase = 3π/16

methods = haskey(args, "methods") ? Symbol.(split(args["methods"], ",")) : methods_run
probs = haskey(args, "p") ? parse.(Float64, split(args["p"], ",")) : probs_run
sample_ids = haskey(args, "samples") ? parse.(Int, split(args["samples"], ",")) : collect(1:nsamples_run)

issubset(methods, (:campspp, :camps)) || error("Only campspp and camps can be replayed, got $methods")
issubset(methods, methods_run) || error("Methods $methods not all in the run ($methods_run)")
issubset(probs, probs_run) || error("p values $probs not all in the run ($probs_run)")
all(in(1:nsamples_run), sample_ids) || error("Sample ids must lie in 1:$nsamples_run")
allunique(methods) && allunique(probs) && allunique(sample_ids) || error("Repeated argument values")

outdir = get(args, "out", joinpath("output", "t11", inject == :none ? runid : "$(runid)_inject-$inject"))
mkpath(outdir)

# == Saved inputs =============================================================
""" Sample id => excited-qubit indices, checked against the bitstring. """
function read_initstates(file, N)
  states = Dict{Int, Vector{Int}}()
  for line in eachline(file)
    (startswith(line, "#") || isempty(strip(line))) && continue
    t = split(line)
    s, bits, inds = parse(Int, t[1]), t[2], parse.(Int, t[3:end])
    length(bits) == N || error("Sample $s: bitstring of length $(length(bits)), expected $N")
    inds == findall(==('1'), collect(bits)) || error("Sample $s: indices don't match the bitstring")
    haskey(states, s) && error("Sample $s listed twice in $file")
    states[s] = inds
  end
  return states
end

""" (p, sample id) => [(gate index, cycle, Pauli string, qubits), ...] """
function read_doped(file)
  doped = Dict{Tuple{Float64, Int}, Vector{Tuple{Int, Int, String, Vector{Int}}}}()
  p = nothing; cur = nothing; expected = 0
  check_count() = isnothing(cur) || length(cur) == expected ||
    error("p = $p: a sample lists $(length(cur)) doped gates, its count line says $expected")
  for line in eachline(file)
    if (m = match(r"^# p = (\S+), μ = ", line)) !== nothing
      check_count(); cur = nothing
      p = parse(Float64, m[1])
    elseif (m = match(r"^# sample (\d+): (\d+) doped gates$", line)) !== nothing
      check_count()
      key = (p, parse(Int, m[1]))
      haskey(doped, key) && error("Doped gates for $key listed twice")
      cur = doped[key] = Tuple{Int, Int, String, Vector{Int}}[]
      expected = parse(Int, m[2])
    elseif !startswith(line, "#") && !isempty(strip(line))
      t = split(line)
      push!(cur, (parse(Int, t[1]), parse(Int, t[2]), String(t[3]), parse.(Int, t[4:end])))
    end
  end
  check_count()
  return doped
end

initstates = read_initstates(initstates_file, N)
doped = read_doped(doped_file)
sort(collect(keys(initstates))) == 1:nsamples_run ||
  error("initstates has samples $(sort(collect(keys(initstates)))), expected 1:$nsamples_run")
for p in probs_run, s in 1:nsamples_run
  haskey(doped, (p, s)) || error("No doped gates saved for p = $p, sample $s")
end

gates, phases_ref = xxz_circuit(ϕ, θ, N/2, N)
layer_ends = layerends(N, N ÷ 2, xxz_circuit)
ntriples = length(phases_ref) ÷ 3

""" The doping mask of one saved sample, after checking each saved gate against the circuit. """
function saved_mask(p, s)
  mask = falses(ntriples)
  list = doped[(p, s)]
  for (k, cycle, str, qinds) in list
    g = gates[k]
    join(string.(g.symbols)) == str && g.qinds == qinds ||
      error("p = $p, sample $s: saved gate $k is $str $qinds, the circuit has $(g.symbols) $(g.qinds)")
    searchsortedfirst(layer_ends, k) == cycle ||
      error("p = $p, sample $s: saved cycle $cycle of gate $k doesn't match layer_ends")
    mask[(k + 2) ÷ 3] = true
  end
  # The doping must be exactly the XX and YY gates of the masked triples
  doped_phases = xy_magic(phases_ref, mask; magicphase = magic_phase)
  findall(doped_phases .!= phases_ref) == first.(list) ||
    error("p = $p, sample $s: saved gates are not whole XX-YY pairs of XXZ triples")
  return mask
end

""" Sensitivity-check defects: `:mask` moves each doped triple one triple later
(the last triple wraps to the first), `:phase` scales the magic phase by 1.05. """
function phases_for(mask)
  m = inject == :mask ? circshift(mask, 1) : mask
  ph = inject == :phase ? 1.05 * magic_phase : magic_phase
  return xy_magic(phases_ref, m; magicphase = ph)
end

# Validate every requested (p, sample) up front, so a bad input fails before any run
masks = Dict((p, s) => saved_mask(p, s) for p in probs, s in sample_ids)

# == Informational: does this Julia regenerate the saved inputs? ==============
# Phase 0 drew the masks from `hash((seed, p))`, which Julia doesn't keep stable
# across versions; the replay never relies on it. If both match, the warm-up sample
# (index 0, not saved) is regenerated the same way; otherwise sample 1 stands in.
regen_states = all(s -> domainwallstate(mt.QubitCAMPS, MersenneTwister(rng_seed + s), N, μ)[2] ==
  initstates[s], 1:nsamples_run)
regen_mask(p, s) = doping_mask(MersenneTwister(hash((rng_seed + s, p))), phases_ref, p)
regen_masks = all(regen_mask(p, s) == saved_mask(p, s) for p in probs_run, s in 1:nsamples_run)
logline("Replay of $runid on $(Threads.nthreads()) thread(s), Julia $VERSION; inject = $inject")
logline("Initial states regenerated from the RNG on this Julia: $(regen_states ? "identical" : "DIFFERENT")")
logline("Doping masks regenerated from hash((seed, p)) on this Julia: $(regen_masks ? "identical" : "DIFFERENT")")
warmup_regen = regen_states && regen_masks
logline("Warm-up sample: " * (warmup_regen ? "index 0, regenerated as in Phase 0" : "sample 1 (index 0 can't be regenerated)"))

# == Outputs ==================================================================
param_info = Dict(hdr)
param_info["methods"] = join(methods, ",")
param_info["magic_probs"] = join(probs, ",")
param_info["sample_ids"] = join(sample_ids, ",")
param_info["replay_of"] = runid
param_info["inject"] = string(inject)
param_info["stack"] = "QubitCAMPS"

campspp_thl_exp = abs(round(Int, log10(campspp_thl)))
o(name) = joinpath(outdir, name)
evs_full = Dict(m => o("evs_prob_tm_$(m)_$(campspp_thl_exp)_full.txt") for m in methods)
bd_full = Dict(m => o("resources_prob_tm_bd_full_$m.txt") for m in methods)
times_full = Dict(m => o("comptime_full_prob_tm_$m.txt") for m in methods)
np_full = o("resources_prob_tm_NP_full_campspp.txt")
pw_full = o("resources_prob_tm_PW_full_campspp.txt")
pc_full = o("resources_prob_tm_PC_full_campspp.txt")
switch_output = o("campspp_switch_prob_tm.txt")
doped_output = o("doped_prob_tm.txt")
initstates_output = o("initstates_prob_tm.txt")
campspp_log = o("comptime_N$(N)_campspp_log.txt")
campssrc_log = o("comptime_N$(N)_campssrc_log.txt")

for m in methods
  initialize_output(evs_full[m], "[layer, expectation value; one block per sample, one set per p]", param_info)
  initialize_output(bd_full[m], "[layer, bond dim.; one block per sample, one set per p]", param_info)
  initialize_output(times_full[m], "[sample, time (s), GC time (s); one block per p]", param_info)
end
if :campspp in methods
  for (f, name) in ((np_full, "n. of Paulis"), (pw_full, "avg. Pauli weight"), (pc_full, "avg. |Pauli coeff.|"))
    initialize_output(f, "[layer, $name; one block per sample, one set per p]", param_info)
  end
  initialize_output(switch_output,
    "[sample, switch gate (last gate applied in CAMPS; PP starts at the next one; $(length(gates)) = never switched), n. of magic qubits k at the switch, PP stop gate ($(length(gates)) = end of circuit); one block per p]",
    param_info)
  open(campspp_log, "w") do _ end
end
:camps in methods && open(campssrc_log, "w") do _ end
initialize_output(doped_output,
  "[doped gates: count line, then one gate position + cycle + Pauli string + qubits per gate; one set per sample, one block per p]",
  param_info)
initialize_output(initstates_output,
  "[sample, initial bitstring (character j is qubit j; 1 = excited), indices of excited qubits]",
  param_info)
open(initstates_output, "a") do io
  for s in sample_ids
    println(io, s, " ", join(j in initstates[s] ? '1' : '0' for j in 1:N), " ", join(initstates[s], " "))
  end
end

# == Main loop ================================================================
init(inds, mask) = (computationalcamps(mt.QubitCAMPS, N, inds), transferredmagnetization(N, inds), phases_for(mask))

for p in probs
  pstr = @sprintf("%.3f", p)
  # The doped gates actually run (after any injected defect), for the input check
  save_doped_gates(doped_output, gates,
    [findall(phases_for(masks[(p, s)]) .!= phases_ref) for s in sample_ids], layer_ends, μ, p)

  if warmup_regen
    inds_wu, mask_wu = domainwallstate(mt.QubitCAMPS, MersenneTwister(rng_seed), N, μ)[2], regen_mask(p, 0)
  else
    inds_wu, mask_wu = initstates[1], saved_mask(p, 1)
  end

  if :campspp in methods
    ψ, tm, phases = init(inds_wu, mask_wu)
    logline("p=$pstr CAMPS-PP warm-up")
    campspp_circuit_dynamics(ψ, campspp_χ, campspp_thl, campspp_Pmax, gates, phases, tm,
      campspp_log; layer_ends = layer_ends)

    evs_s, bds_s, nps_s, pws_s, pcs_s = Vector{Real}[], Vector{Int}[], [], [], []
    times, gctimes, sw_gates, sw_ks, stops = Float64[], Float64[], Int[], Int[], Int[]
    for (j, s) in enumerate(sample_ids)
      ψ, tm, phases = init(initstates[s], masks[(p, s)])
      (evs, tstop, bds, nps, pws, pcs, s_switch, k_switch), time, _, gctime = @timed campspp_circuit_dynamics(
        ψ, campspp_χ, campspp_thl, campspp_Pmax, gates, phases, tm, campspp_log;
        layer_ends = layer_ends, track = true)
      length(evs) < N/2 + 1 && logline("WARNING: CAMPS-PP sample $s for p=$pstr stopped early")
      push!(evs_s, evs); push!(bds_s, bds)
      push!(nps_s, [fill(missing, length(bds) - 1); nps])
      push!(pws_s, [fill(missing, length(bds) - 1); pws])
      push!(pcs_s, [fill(missing, length(bds) - 1); pcs])
      push!(times, time); push!(gctimes, gctime)
      push!(sw_gates, s_switch); push!(sw_ks, k_switch); push!(stops, tstop)
      logline(@sprintf("p=%s CAMPS-PP sample %d (%d/%d, %d%%): %.1f s, switch gate %d, k %d",
        pstr, s, j, length(sample_ids), round(Int, 100j / length(sample_ids)), time, s_switch, k_switch))
    end
    save_full_samples(evs_full[:campspp], μ, p, stack_samples(evs_s))
    save_full_samples(bd_full[:campspp], μ, p, stack_samples(bds_s))
    save_full_samples(np_full, μ, p, stack_samples(nps_s))
    save_full_samples(pw_full, μ, p, stack_samples(pws_s))
    save_full_samples(pc_full, μ, p, stack_samples(pcs_s))
    save_per_sample(times_full[:campspp], μ, p, times, gctimes)
    save_per_sample(switch_output, μ, p, sw_gates, sw_ks, stops)
  end

  if :camps in methods
    ψ, tm, phases = init(inds_wu, mask_wu)
    logline("p=$pstr CAMPS warm-up")
    campssrc_circuit_dynamics(ψ, gates, phases, camps_thl, tm, campssrc_log;
      criterion = camps_crit, strategy = camps_strat, layer_ends = layer_ends)

    evs_s, bds_s, times, gctimes = Vector{Real}[], Vector{Int}[], Float64[], Float64[]
    for (j, s) in enumerate(sample_ids)
      ψ, tm, phases = init(initstates[s], masks[(p, s)])
      (_, _, evs, bds), time, _, gctime = @timed campssrc_circuit_dynamics(
        ψ, gates, phases, camps_thl, tm, campssrc_log;
        criterion = camps_crit, strategy = camps_strat, layer_ends = layer_ends, track = true)
      push!(evs_s, evs); push!(bds_s, bds); push!(times, time); push!(gctimes, gctime)
      logline(@sprintf("p=%s CAMPS sample %d (%d/%d, %d%%): %.1f s, max bond dim. %d",
        pstr, s, j, length(sample_ids), round(Int, 100j / length(sample_ids)), time, maximum(bds)))
    end
    save_full_samples(evs_full[:camps], μ, p, stack_samples(evs_s))
    save_full_samples(bd_full[:camps], μ, p, stack_samples(bds_s))
    save_per_sample(times_full[:camps], μ, p, times, gctimes)
  end
end
logline("Done; outputs in $outdir")
