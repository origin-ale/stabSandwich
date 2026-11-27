using CampsPP
using CampsPP.DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp

using Statistics
using Random: MersenneTwister

@isdefined(ScriptSetup) || include(joinpath(@__DIR__, "scripts", "ScriptSetup.jl"))
using .ScriptSetup

# Method comparison on the self-dual kicked Ising model (KIM) near its Clifford point, with
# periodic boundaries, from the Néel state (test/CONVENTIONS.md §2 "Kicked Ising"). A run is
# one (sample, site j): ⟨Z_j⟩ at t = 0 and after every layer, until the method reaches its
# budget (checked at layer ends) or the cap t_max.

# == Parameters ===============================================================
# Override from the command line as key=value (`help` lists them); lists are
# comma-separated, without spaces, e.g.
#   julia --project=. Comptime_kim.jl methods=campspp,pp εs=0,0.05 sites=1,2
params = @parameters begin
  # Methods to run, in output-column order.
  # Available: :campspp (CAMPS-PP), :camps (CAMPS), :mps (MPS), :pp (PP)
  methods = [:campspp, :pp]

  N = 46 # System size (periodic chain)
  # Layer cap: a run that never reaches its budget stops after t_max layers. Needed
  # because at ε = 0 PP keeps one Pauli string and CAMPS-PP never switches.
  t_max = 2N

  # Deviations ε from the Clifford point of the Z field, strictly increasing
  εs = [0.0, 0.01, 0.02, 0.05, 0.1, 0.2]
  # Sites j whose ⟨Z_j⟩ is recorded, one run each
  sites = collect(1:N)
  # :uniform: ε_j = ε on every qubit (needs samples = 1). :random: sample i draws, at each
  # ε, u ∈ [−1, 1]^N from seed hash((rng_seed + i, ε)) and uses ε_j = ε·u_j (u independent
  # across ε, shared across methods and sites; sample 0 is the untimed warm-up)
  eps_mode = :uniform
  samples = 1
  rng_seed = 100

  # CAMPS-PP (sandwich, ie. s-)
  campspp_χ = 64 # s-CAMPS bond dimension that triggers s-PP
  campspp_thl = 1e-10 # Threshold for truncation during s-PP (> 0)
  campspp_Pmax = 100_000 # Budget: stop once s-PP has ≥ this many Pauli strings

  # CAMPS
  camps_thl = 1e-12 # CAMPS SVD truncation threshold
  camps_crit = :chi3 # CAMPS entanglement criterion (:chi3 or :entangle)
  camps_strat = :snake # CAMPS disentangler strategy (:snake or :brickwork)
  camps_χmax = 256 # Budget: stop once the bond dimension is ≥ this

  # MPS
  mps_thl = 1e-7 # MPS SVD truncation threshold
  mps_χmax = 256 # Budget: stop once the bond dimension is ≥ this

  # Pauli propagation
  pp_thl = 1e-10 # Threshold for truncation during PP (> 0)
  pp_Pmax = 100_000 # Budget: stop once PP has ≥ this many Pauli strings

  # BLAS and Strided threads; runs go one after another (see Comptime_prob_tm.jl)
  blas_threads = 1

  outdir = "output/Comptime_kim"
end

# == Validation (before any output) ===========================================

all_methods = [:campspp, :camps, :mps, :pp]
let unknown = setdiff(methods, all_methods)
  isempty(unknown) || error("Unrecognised method(s): $unknown")
end
isempty(methods) && error("No methods selected to run")
allunique(methods) || error("Repeated method(s) in $methods")
eps_mode in (:uniform, :random) ||
  error("eps_mode must be :uniform or :random, got :$eps_mode")
samples >= 1 || error("Need at least 1 sample, got samples = $samples")
eps_mode == :uniform && samples != 1 &&
  error("eps_mode = :uniform needs samples = 1 (every sample would be the same), got samples = $samples")
isempty(εs) && error("No εs given")
issorted(εs; lt = <=) || error("εs must be strictly increasing: $εs")
isempty(sites) && error("No sites given")
allunique(sites) || error("Repeated site(s) in sites = $sites")
all(in(1:N), sites) || error("sites must lie in 1:N = 1:$N, got $sites")
:campspp in methods && !(campspp_thl > 0) &&
  error("campspp_thl must be > 0 (PP takes the Clifford gates as generic rotations), got $campspp_thl")
:pp in methods && !(pp_thl > 0) &&
  error("pp_thl must be > 0 (PP takes the Clifford gates as generic rotations), got $pp_thl")
t_max >= 1 || error("t_max must be ≥ 1, got $t_max")

budget_names = Dict(:campspp => "campspp_Pmax", :camps => "camps_χmax",
  :mps => "mps_χmax", :pp => "pp_Pmax")
budgets = Dict(:campspp => campspp_Pmax, :camps => camps_χmax,
  :mps => mps_χmax, :pp => pp_Pmax)
for m in methods
  b = budgets[m]
  (b isa Integer && !(b isa Bool) && b >= 2) ||
    error("$(budget_names[m]) must be an integer ≥ 2, got $b")
end

εs = Float64.(εs)
setup_threads(blas_threads)

# Threshold exponent of each method, appended to its expectation-value file names ("exact"
# for a threshold ≤ 0, which truncates nothing)
thl_exp(thl) = thl <= 0 ? "exact" : abs(round(Int, log10(thl)))
thl_methods = Dict(:campspp => campspp_thl, :camps => camps_thl, :mps => mps_thl, :pp => pp_thl)

param_info = param_header(params; exclude = (:outdir,))

# == Output files =============================================================

mkpath(outdir)
output = joinpath(outdir, "comptime_kim.txt")
gctime_output = joinpath(outdir, "gctime_kim.txt")
layers_output = joinpath(outdir, "layers_kim.txt")
runs_output = joinpath(outdir, "runs_kim.txt")
eps_output = joinpath(outdir, "eps_kim.txt")
switch_output = joinpath(outdir, "campspp_switch_kim.txt")
resources_prefix = joinpath(outdir, "resources_kim_")
evs_prefix = joinpath(outdir, "evs_kim_")

logs = Dict(m => joinpath(outdir, "kim_N$(N)_$(m == :camps ? "campssrc" : m)_log.txt")
  for m in all_methods)

# Runs, sample-major: run r = (i − 1)·|sites| + (position of j in sites)
R = samples * length(sites)
run_sample = [i for i in 1:samples for _ in sites]
run_site = [j for _ in 1:samples for j in sites]
initialize_output(runs_output, "[run, sample, site j]", param_info)
save_columns(runs_output, collect(1:R), run_sample, run_site)

eps_mode == :random && initialize_output(eps_output,
  "[sample i, u_1 … u_N (ε_j = ε·u_j); one block per ε]", param_info)

# Per-method summary columns, one entry per ε: mean, std. err., q1, median, q3, min, max
stat_names = (:mean, :err, :q1, :med, :q3, :min, :max)
new_stats() = Dict(m => Dict(s => Float64[] for s in stat_names) for m in methods)
time_stats, gc_stats, layer_stats = new_stats(), new_stats(), new_stats()

"Push the seven statistics of `xs` (every run) onto `stats[m]`."
function push_stats!(stats, m, xs)
  q1, med, q3 = quantiles(xs)
  vals = (mean(xs), std(xs) / sqrt(length(xs)), q1, med, q3, minimum(xs), maximum(xs))
  for (s, v) in zip(stat_names, vals)
    push!(stats[m][s], v)
  end
end

"Rewrite a summary table: x column ε (the εs done so far), seven columns per method."
function write_summary(file, what, stats, εs_done)
  initialize_output(file,
    "[ε, then per method: $what mean, std. err., q1, median, q3, min, max; over every run]",
    param_info)
  save_columns(file, εs_done,
    (stats[m][s] for m in methods for s in stat_names)...)
end

# Resource files, as in Comptime_prob_tm.jl
bd_methods = [m for m in methods if m in (:campspp, :camps)]
pp_methods = [m for m in methods if m in (:campspp, :pp)]
resource_kinds = [("bd", "bond dim.", bd_methods), ("NP", "n. of Paulis", pp_methods),
  ("PW", "avg. Pauli weight", pp_methods), ("PC", "avg. |Pauli coeff.|", pp_methods)]
res_outputs = Dict((k, m) => resources_prefix * "$(k)_$m.txt"
  for (k, _, ms) in resource_kinds for m in ms)
res_full_outputs = Dict((k, m) => resources_prefix * "$(k)_full_$m.txt"
  for (k, _, ms) in resource_kinds for m in ms)
res_peak_outputs = Dict((k, m) => resources_prefix * "$(k)_peak_$m.txt"
  for (k, _, ms) in resource_kinds for m in ms)
for (k, name, ms) in resource_kinds, m in ms
  initialize_output(res_outputs[(k, m)], "[layer, mean $name, std. err., max-run $name, n. of runs, q1, median, q3, min, max; one block per ε]", param_info)
  initialize_output(res_full_outputs[(k, m)], "[layer, $name; one block per run (# sample r), one set per ε]", param_info)
  initialize_output(res_peak_outputs[(k, m)], "[run, peak $name; one block per ε]", param_info)
end

times_full_outputs = Dict(m => joinpath(outdir, "comptime_full_kim_$m.txt") for m in methods)
layers_full_outputs = Dict(m => joinpath(outdir, "layers_kim_$m.txt") for m in methods)
evs_outputs = Dict(m => evs_prefix * "$(m)_$(thl_exp(thl_methods[m])).txt" for m in methods)
evs_full_outputs = Dict(m => evs_prefix * "$(m)_$(thl_exp(thl_methods[m]))_full.txt" for m in methods)
for m in methods
  initialize_output(times_full_outputs[m], "[run, time (s), GC time (s); one block per ε]", param_info)
  initialize_output(layers_full_outputs[m], "[run, stop layer l, s (1 = budget $(budget_names[m]) reached at l, 0 = stopped on the cap t_max); one block per ε]", param_info)
  initialize_output(evs_outputs[m], "[layer, then per site in sites order ($(join(sites, ","))): mean ⟨Z_j⟩, std. err., over the runs that recorded the layer; one block per ε]", param_info)
  initialize_output(evs_full_outputs[m], "[layer, ⟨Z_j⟩; one block per run (# sample r), one set per ε]", param_info)
end

n_gates = last(kicked_ising_layerends(N, t_max))
if :campspp in methods
  initialize_output(switch_output,
    "[run, switch gate (last gate applied in CAMPS; PP starts at the next one; $n_gates = never switched), n. of magic qubits k at the switch, PP stop gate; one block per ε]",
    param_info)
end

# == Helpers ==================================================================

"u ∈ [−1, 1]^N of sample i at ε (`:random`; i = 0 is the warm-up)."
sample_u(i, ε) = 2 .* rand(MersenneTwister(hash((rng_seed + i, ε))), N) .- 1

"ε_j of sample i at ε."
sample_eps(i, ε) = eps_mode == :uniform ? fill(ε, N) : ε .* sample_u(i, ε)

neel_ones = collect(2:2:N)

function neel_mps()
  strs = fill("0", N)
  strs[neel_ones] .= "1"
  return MPS(siteinds("Qubit", N), strs)
end

"""
Profile block: per layer 0 … (deepest recorded), for each site in `sites` order, the mean of
⟨Z_j⟩ and its standard error over the runs at site j that recorded the layer (NaN mean if
none, NaN error below two).
"""
function save_profile(file, label, evs_runs)
  deepest = maximum(length, evs_runs) - 1
  open(file, "a") do io
    println(io, CampsPP.label_header(label, nothing))
    for l in 0:deepest
      fields = Any[l]
      for j in sites
        xs = Float64[evs_runs[r][l + 1] for r in 1:R
                     if run_site[r] == j && length(evs_runs[r]) > l]
        push!(fields, isempty(xs) ? NaN : mean(xs))
        push!(fields, length(xs) < 2 ? NaN : std(xs) / sqrt(length(xs)))
      end
      println(io, join(string.(fields), " "))
    end
    print(io, "\n\n")
  end
end

""" Indices of the runs whose ⟨Z_j⟩ in `as` and `bs` differ by more than `atol` at some \
layer that both recorded. """
mismatched_runs(as, bs; atol) = [r for r in eachindex(as, bs)
  if !all(isapprox.(as[r][1:min(end, length(bs[r]))], bs[r][1:min(end, length(as[r]))]; atol = atol))]

# == Main loop ================================================================
layer_ends = kicked_ising_layerends(N, t_max)

for (iε, ε) in enumerate(εs)
  label = "ε" => ε
  evs_by_method = Dict{Symbol, Vector{Vector{Float64}}}()

  eps_runs = [sample_eps(i, ε) for i in 1:samples]
  eps_wu = sample_eps(0, ε)
  if eps_mode == :random
    open(eps_output, "a") do io
      println(io, CampsPP.label_header(label, nothing))
      for i in 1:samples
        println(io, i, " ", join(string.(sample_u(i, ε)), " "))
      end
      print(io, "\n\n")
    end
  end
  circuit(i) = kicked_ising_circuit(i == 0 ? eps_wu : eps_runs[i], t_max, N)

  for m in methods
    prog = ProgressLog("ε=$ε $m run", R; frac = 0)
    times, gctimes = Float64[], Float64[]
    stop_layers, stop_flags = Int[], Int[]
    evs_runs = Vector{Float64}[]
    res = Dict(k => AbstractVector[] for (k, _, _) in resource_kinds)
    switch_gates, switch_ks, stop_gates = Int[], Int[], Int[]
    B = budgets[m]

    # One run of method m for (sample i, site j): the driver's result, timed with @timed.
    # The circuit, observable and initial state are built outside the timing (as in
    # Comptime_prob_tm.jl).
    function run_method(i, j)
      gates, phases = circuit(i)
      obs = local_z(N, j)
      if m == :campspp
        ψ = computationalcamps(N, neel_ones)
        f = () -> campspp_circuit_dynamics(ψ, campspp_χ,
          campspp_thl, B, gates, phases, obs, logs[m]; layer_ends = layer_ends, track = true)
      elseif m == :camps
        ψ = computationalcamps(N, neel_ones)
        f = () -> campssrc_circuit_dynamics(ψ, gates, phases,
          camps_thl, obs, logs[m]; criterion = camps_crit, strategy = camps_strat,
          layer_ends = layer_ends, track = true, χmax = B)
      elseif m == :mps
        ψ = neel_mps()
        f = () -> mps_circuit_dynamics(ψ, gates, phases, mps_thl, obs, logs[m];
          layer_ends = layer_ends, χmax = B)
      else
        ones_pp = copy(neel_ones)
        f = () -> pauliprop_circuit_dynamics(ones_pp, gates, phases, pp_thl, B, obs,
          logs[m]; layer_ends = layer_ends, track = true)
      end
      return @timed f()
    end

    run_method(0, first(sites))    # untimed warm-up

    for r in 1:R
      out, time, _, gctime = run_method(run_sample[r], run_site[r])
      if m == :campspp
        evs, tstop, bds, nps, pws, pcs, s_switch, k_switch = out
        # PP starts one past the last CAMPS-recorded layer; pad so that row k is layer k
        pad = fill(missing, length(bds) - 1)
        push!(res["bd"], bds)
        push!(res["NP"], [pad; nps])
        push!(res["PW"], [pad; pws])
        push!(res["PC"], [pad; pcs])
        push!(switch_gates, s_switch); push!(switch_ks, k_switch); push!(stop_gates, tstop)
        # The hand-over entry (1 string) is < B, so a run that never switched has s = 0
        reached = length(nps) > 1 && last(nps) >= B
      elseif m == :camps
        _, _, evs, bds = out
        push!(res["bd"], bds)
        reached = last(bds) >= B
      elseif m == :mps
        ψ, evs = out
        reached = maxlinkdim(ψ) >= B
      else
        _, evs, nps, pws, pcs = out
        push!(res["NP"], nps); push!(res["PW"], pws); push!(res["PC"], pcs)
        reached = last(nps) >= B
      end
      push!(evs_runs, Float64.(evs))
      push!(times, time); push!(gctimes, gctime)
      push!(stop_layers, length(evs) - 1)
      push!(stop_flags, reached ? 1 : 0)
      log_progress!(prog, r)
    end

    ncap = count(==(0), stop_flags)
    log_message("ε=$ε $m: $ncap of $R runs stopped on the cap t_max = $t_max without reaching $(budget_names[m]) = $B")

    push_stats!(time_stats, m, times)
    push_stats!(gc_stats, m, gctimes)
    push_stats!(layer_stats, m, Float64.(stop_layers))
    save_per_sample(times_full_outputs[m], nothing, label, times, gctimes)
    save_per_sample(layers_full_outputs[m], nothing, label, stop_layers, stop_flags)
    save_profile(evs_outputs[m], label, evs_runs)
    save_full_samples(evs_full_outputs[m], nothing, label, stack_samples(evs_runs))
    for (k, _, ms) in resource_kinds
      m in ms || continue
      save_stats_maxcol(res_outputs[(k, m)], res[k], nothing, label)
      save_full_samples(res_full_outputs[(k, m)], nothing, label, stack_samples(res[k]))
      save_peaks(res_peak_outputs[(k, m)], res[k], nothing, label)
    end
    m == :campspp &&
      save_per_sample(switch_output, nothing, label, switch_gates, switch_ks, stop_gates)
    evs_by_method[m] = evs_runs
  end

  # == Cross-checks: CAMPS-PP against each other method, on the layers both recorded
  if haskey(evs_by_method, :campspp)
    for (m, name) in ((:camps, "CAMPS"), (:mps, "MPS"), (:pp, "PP"))
      haskey(evs_by_method, m) || continue
      a, b = evs_by_method[:campspp], evs_by_method[m]
      bad = mismatched_runs(a, b; atol = 1e-6)
      isempty(bad) && continue
      log_message("WARNING: ε = $ε CAMPS-PP and $name do not match (runs $bad)")
      println("-"^32)
      for r in bad
        println("$(a[r])\n$(b[r])\n")
      end
      println("-"^32)
    end
  end

  write_summary(output, "time (s)", time_stats, εs[1:iε])
  write_summary(gctime_output, "GC time (s)", gc_stats, εs[1:iε])
  write_summary(layers_output, "stop layer", layer_stats, εs[1:iε])
end
