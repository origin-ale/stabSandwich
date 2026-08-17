using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps
using QuantumClifford

using Statistics
using Random: MersenneTwister
using ProgressMeter
using Strided
using LinearAlgebra

Strided.disable_threads()
nthr = Threads.nthreads()

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

# All three methods are run in their *exact* settings, so any disagreement between
# them is a bug and not a truncation artefact:
#   CAMPS — analytical disentangling, intrinsically exact
#   MPS   — cutoff 0 and naive (pairwise) MPO·MPS contraction
#   PP    — min_abs_coeff 0, every other truncation left at Inf

function circuit_camps(gates, phases, P)
  ψ, _ = evolve(cmps.CAMPS(P.nqubits), length(phases), gates, phases)
  return real(cmps.expectation(ψ, P))
end

function circuit_mps(gates, phases, P; cutoff = 0, alg = "naive")
  N = P.nqubits
  I = PauliOperator(0x0, fill(false,N), fill(false,N))

  gate_paulis = cmps.PauliOperator.(getpauli.(gates, N))

  sites = siteinds("Qubit", N)
  states = ["Up" for _ in sites]
  ψ = MPS(sites, states)

  for i in eachindex(gate_paulis)
    p = gate_paulis[i]
    ϕ = phases[i]
    g = cmps.PauliSum([cos(ϕ), im * sin(ϕ)], Stabilizer([I,p]))
    ψ = apply(ψ, g; cutoff = cutoff, alg = alg)
  end

  P_cmps = cmps.PauliSum(P)
  return real(ITensorMPS.expect(ψ, P_cmps))
end

function circuit_pp(gates, phases, P, thl)
  propagated = pp.propagate(gates, P, -2 .* phases; min_abs_coeff = thl)
  return real(pp.overlapwithzero(propagated))
end

# The pp.PauliRotation -> QuantumClifford.PauliOperator conversion that both the CAMPS
# and the MPS path perform internally, isolated so that it can be timed on its own.
convert_gates(gates, N) = PauliOperator.(getpauli.(gates, N))

# == Parameters ===============================================================
# Methods to run, in output-column order. Comment out any you want to skip.
# Available: :camps (CAMPS), :pp (PP), :mps (MPS)
methods = [:camps, :pp, :mps]

Nsamples = 100
Nrange = [2, 4, 6, 8, 10, 12, 16, 22, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 44, 64, 90, 128]
Nmax_mps_run = 25 # Exact MPS is only run for N ≤ Nmax_mps_run
Nmax_pp_run = 70 # Exact PP is only run for N ≤ Nmax_pp_run

# Base RNG seed: sample i uses seed rng_seed + i
rng_seed = 100
# Warm-up circuits run (untimed) before each N, on seeds rng_seed, rng_seed-1, …
# More than one because CAMPS branches on the Clifford nature of each rotation,
# so a single circuit does not reach every code path
warmup_samples = 10

# MPS
cutoff_mps = 0.0 # Exact: no SVD truncation
alg_mps = "naive" # Exact: pairwise MPO·MPS contraction

# Pauli propagation
thl_pp = 0.0 # Exact: no coefficient truncation

# Agreement tolerance of the cross-method check
atol_check = 1e-6

# =============================================================================

all_methods = [:camps, :mps, :pp]
let unknown = setdiff(methods, all_methods)
  isempty(unknown) || error("Unrecognised method(s): $unknown")
end
isempty(methods) && error("No methods selected to run")

# Per-method N ranges; both are prefixes of Nrange, so the padding of short
# columns in `save_columns` lines the output up
Nrange_mps = filter(≤(Nmax_mps_run), Nrange)
Nrange_pp = filter(≤(Nmax_pp_run), Nrange)
# The gate conversion is timed at every N, like CAMPS
Nranges = Dict(:camps => Nrange, :mps => Nrange_mps, :pp => Nrange_pp,
  :convert => Nrange)

obs_string = "Z₁"
obsname = "⟨$obs_string⟩ computation time (s)"

out_full = "output/comptimes_3_squarecirc_full.txt"
out_avgs = "output/comptimes_3_squarecirc_avgs.txt"
out_diffs_full = "output/timediffs_3_squarecirc_full.txt"
out_diffs_avgs = "output/timediffs_3_squarecirc_avgs.txt"

param_info = Dict(
  "methods" => methods,
  "Nrange" => Nrange,
  "Nsamples" => Nsamples,
  "Nmax_mps_run" => Nmax_mps_run,
  "Nmax_pp_run" => Nmax_pp_run,
  "rng_seed" => rng_seed,
  "warmup_samples" => warmup_samples,
  "cutoff_mps" => cutoff_mps,
  "alg_mps" => alg_mps,
  "thl_pp" => thl_pp,
  "obs" => obs_string)

# Timings, one entry per N; :convert is the gate-conversion sub-timing, run at every N
times = Dict(m => Vector{Float64}[] for m in [methods; :convert])
# Per-realization CAMPS - PP wall-clock time, one entry per N where both methods ran
timediffs = Vector{Float64}[]
Nrange_diffs = Int[]
# Runs whose warm-up failed to keep the compiler out of the timing, per method
warmup_misses = Dict{Symbol, Int}()

""" Warn if compilation landed inside a timed run, ie. if the warm-up missed. """
function check_warm(st, m, N, i)
  iszero(st.compile_time) && return
  printstyled("WARNING: N = $N, sample $i, $m: compile_time = $(st.compile_time) s \
($(round(100 * st.compile_time / st.time; digits = 1))% of $(st.time) s), warm-up missed\n";
    color = :yellow)
  warmup_misses[m] = get(warmup_misses, m, 0) + 1
end

""" Warn if the methods that ran on this realization disagree. All of them are
exact, so any disagreement beyond `atol_check` is a bug. """
function check_evs(evs, N, i)
  ms = [m for m in methods if haskey(evs, m)]
  length(ms) < 2 && return
  ref = first(ms)
  for m in ms[2:end]
    if !isapprox(evs[ref], evs[m]; atol = atol_check)
      printstyled("WARNING: N = $N, sample $i: $ref and $m do not match \
($(evs[ref]) vs $(evs[m]), difference $(evs[ref] - evs[m]))\n"; color = :yellow)
    end
  end
end

initialize_output(out_full,
  "$obsname, $Nsamples samples [N, then one column per sample; one block per \
method in the order $(join(methods, ", ")), then one for the gate conversion]",
  param_info)
initialize_output(out_avgs,
  "$obsname averages [N, then per method: mean time (s), std. err., in the order \
$(join(methods, ", ")); last two columns are the gate conversion]",
  param_info)
initialize_output(out_diffs_full,
  "CAMPS - PP computation time (s), $Nsamples samples [N, then one column per sample]",
  param_info)
initialize_output(out_diffs_avgs,
  "CAMPS - PP computation time (s) [N, mean, std. err.]",
  param_info)

printstyled("Getting exact $(join(methods, ", ")) computation times of ⟨$obs_string⟩ for \
square circuits, N = $Nrange.\n\
MPS exact for N ≤ $Nmax_mps_run (cutoff $cutoff_mps, alg $alg_mps), \
PP exact for N ≤ $Nmax_pp_run (min_abs_coeff $thl_pp).\n\
Nsamples = $Nsamples, $nthr threads.\n"; color = :cyan)
prog = Progress(length(Nrange) * Nsamples; desc = "Computing…")

for N in Nrange
  P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
  active = [m for m in methods if N in Nranges[m]]

  # Warm-up on circuits disjoint from the sampled ones, results discarded:
  # compilation must not land in the timings of sample 1
  for w in 0:warmup_samples-1
    gates_wu, phases_wu = rotation_circuit(MersenneTwister(rng_seed - w), N, N)
    convert_gates(gates_wu, N)
    :camps in active && circuit_camps(gates_wu, phases_wu, P)
    :mps in active && circuit_mps(gates_wu, phases_wu, P;
      cutoff = cutoff_mps, alg = alg_mps)
    :pp in active && circuit_pp(gates_wu, phases_wu, P, thl_pp)
  end

  times_N = Dict(m => Float64[] for m in [active; :convert])

  for i in 1:Nsamples
    # One circuit per (N, sample), shared by every method: the comparison is on
    # identical realizations by construction
    gates, phases = rotation_circuit(MersenneTwister(rng_seed + i), N, N)
    evs = Dict{Symbol, Float64}()

    st = @timed convert_gates(gates, N)
    check_warm(st, :convert, N, i)
    push!(times_N[:convert], st.time)

    if :camps in active
      st = @timed circuit_camps(gates, phases, P)
      check_warm(st, :camps, N, i)
      push!(times_N[:camps], st.time)
      evs[:camps] = st.value
    end
    if :mps in active
      st = @timed circuit_mps(gates, phases, P; cutoff = cutoff_mps, alg = alg_mps)
      check_warm(st, :mps, N, i)
      push!(times_N[:mps], st.time)
      evs[:mps] = st.value
    end
    if :pp in active
      st = @timed circuit_pp(gates, phases, P, thl_pp)
      check_warm(st, :pp, N, i)
      push!(times_N[:pp], st.time)
      evs[:pp] = st.value
    end

    check_evs(evs, N, i)

    next!(prog, showvalues = [("N", N), ("sample", i)])
  end

  for m in keys(times_N)
    push!(times[m], times_N[m])
  end
  # Paired per-realization time difference: both methods ran the same circuit,
  # so the difference is taken sample by sample rather than between the averages
  if :camps in active && :pp in active
    push!(timediffs, times_N[:camps] .- times_N[:pp])
    push!(Nrange_diffs, N)
  end
end

for m in [methods; :convert]
  isempty(times[m]) && continue
  save_rows(out_full, Nranges[m], times[m]; blockend = true)
end

stderr_of(xs) = std(xs) / sqrt(length(xs))

""" Spread a method's per-N values over the whole `Nrange`, writing NaN where
that method did not run. `save_columns` would otherwise pad a short column with
*empty* fields, and gnuplot collapses consecutive whitespace: every column after
a method with a short range would shift left on those rows. """
function over_Nrange(m, vals)
  byN = Dict(zip(Nranges[m], vals))
  return [get(byN, N, NaN) for N in Nrange]
end

avg_cols = AbstractVector[]
for m in [methods; :convert]
  isempty(times[m]) && continue
  push!(avg_cols, over_Nrange(m, mean.(times[m])))
  push!(avg_cols, over_Nrange(m, stderr_of.(times[m])))
end
save_columns(out_avgs, Nrange, avg_cols...)

if !isempty(timediffs)
  save_rows(out_diffs_full, Nrange_diffs, timediffs)
  save_columns(out_diffs_avgs,
               Nrange_diffs,
               mean.(timediffs),
               stderr_of.(timediffs))
end

if isempty(warmup_misses)
  printstyled("Warm-up clean: no compilation inside any timed run.\n"; color = :green)
else
  printstyled("WARNING: compilation landed inside timed runs: \
$(join(("$m: $n" for (m, n) in warmup_misses), ", "))\n"; color = :yellow)
end
