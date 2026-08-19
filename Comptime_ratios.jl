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

# Both methods are run in their *exact* settings, so any disagreement between
# them is a bug and not a truncation artefact:
#   CAMPS — analytical disentangling, intrinsically exact
#   PP    — min_abs_coeff 0, every other truncation left at Inf

function circuit_camps(gates, phases, P)
  ψ, _ = evolve(cmps.CAMPS(P.nqubits), length(phases), gates, phases)
  return real(cmps.expectation(ψ, P))
end

function circuit_pp(gates, phases, P, thl)
  propagated = pp.propagate(gates, P, -2 .* phases; min_abs_coeff = thl)
  return real(pp.overlapwithzero(propagated))
end

# The pp.PauliRotation -> QuantumClifford.PauliOperator conversion that the CAMPS
# path performs internally, isolated so that it can be timed on its own.
convert_gates(gates, N) = PauliOperator.(getpauli.(gates, N))

# == Parameters ===============================================================
# Methods to run, in output-column order. Comment out any you want to skip.
# Available: :camps (CAMPS), :pp (PP)
methods = [:camps, :pp]

N = 40
Nsamples = 100

# Aspect ratios t/N to sweep, at fixed N
ratiostart = 0.7
ratioend = 1.2
ratiostep = 0.05

# Base RNG seed: sample i uses seed rng_seed + i
rng_seed = 100
# Warm-up circuits run (untimed) before each ratio, on seeds rng_seed, rng_seed-1, …
# More than one because CAMPS branches on the Clifford nature of each rotation,
# so a single circuit does not reach every code path
warmup_samples = 10

# Pauli propagation
thl_pp = 0.0 # Exact: no coefficient truncation

# Agreement tolerance of the cross-method check
atol_check = 1e-6

# =============================================================================

all_methods = [:camps, :pp]
let unknown = setdiff(methods, all_methods)
  isempty(unknown) || error("Unrecognised method(s): $unknown")
end
isempty(methods) && error("No methods selected to run")

ratios = collect(ratiostart:ratiostep:ratioend)
Nratios = length(ratios)
# Rotation counts of each swept ratio; rounding can collapse two ratios onto the
# same t, so the mapping is reported in the output header
ts = [Int(round(ratio * N)) for ratio in ratios]

obs_string = "Z₁"
obsname = "⟨$obs_string⟩ computation time (s)"

out_full = "output/comptimes_ratios_$(N)_full.txt"
out_avgs = "output/comptimes_ratios_$(N)_avgs.txt"
out_diffs_full = "output/timediffs_ratios_$(N)_full.txt"
out_diffs_avgs = "output/timediffs_ratios_$(N)_avgs.txt"

param_info = Dict(
  "methods" => methods,
  "N" => N,
  "ratios" => ratios,
  "ts" => ts,
  "Nsamples" => Nsamples,
  "rng_seed" => rng_seed,
  "warmup_samples" => warmup_samples,
  "thl_pp" => thl_pp,
  "obs" => obs_string)

# Timings, one entry per ratio; :convert is the gate-conversion sub-timing
times = Dict(m => Vector{Float64}[] for m in [methods; :convert])
# Per-realization CAMPS - PP wall-clock time, one entry per ratio
timediffs = Vector{Float64}[]
# Runs whose warm-up failed to keep the compiler out of the timing, per method
warmup_misses = Dict{Symbol, Int}()

""" Warn if compilation landed inside a timed run, ie. if the warm-up missed. """
function check_warm(st, m, ratio, i)
  iszero(st.compile_time) && return
  printstyled("WARNING: t/N = $ratio, sample $i, $m: compile_time = $(st.compile_time) s \
($(round(100 * st.compile_time / st.time; digits = 1))% of $(st.time) s), warm-up missed\n";
    color = :yellow)
  warmup_misses[m] = get(warmup_misses, m, 0) + 1
end

""" Warn if the methods that ran on this realization disagree. All of them are
exact, so any disagreement beyond `atol_check` is a bug. """
function check_evs(evs, ratio, i)
  ms = [m for m in methods if haskey(evs, m)]
  length(ms) < 2 && return
  ref = first(ms)
  for m in ms[2:end]
    if !isapprox(evs[ref], evs[m]; atol = atol_check)
      printstyled("WARNING: t/N = $ratio, sample $i: $ref and $m do not match \
($(evs[ref]) vs $(evs[m]), difference $(evs[ref] - evs[m]))\n"; color = :yellow)
    end
  end
end

initialize_output(out_full,
  "$obsname, $Nsamples samples [t/N, then one column per sample; one block per \
method in the order $(join(methods, ", ")), then one for the gate conversion]",
  param_info)
initialize_output(out_avgs,
  "$obsname averages [t/N, then per method: mean time (s), std. err., in the order \
$(join(methods, ", ")); last two columns are the gate conversion]",
  param_info)
initialize_output(out_diffs_full,
  "CAMPS - PP computation time (s), $Nsamples samples [t/N, then one column per sample]",
  param_info)
initialize_output(out_diffs_avgs,
  "CAMPS - PP computation time (s) [t/N, mean, std. err.]",
  param_info)

printstyled("Getting exact $(join(methods, ", ")) computation times of ⟨$obs_string⟩ for \
N = $N t-rotation circuits, $Nratios aspect ratios t/N = $ratios (t = $ts).\n\
PP exact (min_abs_coeff $thl_pp).\n\
Nsamples = $Nsamples, $nthr threads.\n"; color = :cyan)
prog = Progress(Nratios * Nsamples; desc = "Computing…")

# N is fixed and every method runs at every ratio, so both are loop-invariant
P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
active = methods

for (r, ratio) in enumerate(ratios)
  t = ts[r]

  # Warm-up on circuits disjoint from the sampled ones, results discarded:
  # compilation must not land in the timings of sample 1. Repeated at every
  # ratio because t changes the circuit length, and so which code paths are hit
  for w in 0:warmup_samples-1
    gates_wu, phases_wu = rotation_circuit(MersenneTwister(rng_seed - w), t, N)
    convert_gates(gates_wu, N)
    :camps in active && circuit_camps(gates_wu, phases_wu, P)
    :pp in active && circuit_pp(gates_wu, phases_wu, P, thl_pp)
  end

  times_r = Dict(m => Float64[] for m in [active; :convert])

  for i in 1:Nsamples
    # One circuit per (ratio, sample), shared by every method: the comparison is
    # on identical realizations by construction
    gates, phases = rotation_circuit(MersenneTwister(rng_seed + i), t, N)
    evs = Dict{Symbol, Float64}()

    st = @timed convert_gates(gates, N)
    check_warm(st, :convert, ratio, i)
    push!(times_r[:convert], st.time)

    if :camps in active
      st = @timed circuit_camps(gates, phases, P)
      check_warm(st, :camps, ratio, i)
      push!(times_r[:camps], st.time)
      evs[:camps] = st.value
    end
    if :pp in active
      st = @timed circuit_pp(gates, phases, P, thl_pp)
      check_warm(st, :pp, ratio, i)
      push!(times_r[:pp], st.time)
      evs[:pp] = st.value
    end

    check_evs(evs, ratio, i)

    next!(prog, showvalues = [("t/N", ratio), ("t", t), ("sample", i)])
  end

  for m in keys(times_r)
    push!(times[m], times_r[m])
  end
  # Paired per-realization time difference: both methods ran the same circuit,
  # so the difference is taken sample by sample rather than between the averages
  if :camps in active && :pp in active
    push!(timediffs, times_r[:camps] .- times_r[:pp])
  end
end

for m in [methods; :convert]
  isempty(times[m]) && continue
  save_rows(out_full, ratios, times[m]; blockend = true)
end

stderr_of(xs) = std(xs) / sqrt(length(xs))

avg_cols = AbstractVector[]
for m in [methods; :convert]
  isempty(times[m]) && continue
  push!(avg_cols, mean.(times[m]))
  push!(avg_cols, stderr_of.(times[m]))
end
save_columns(out_avgs, ratios, avg_cols...)

if !isempty(timediffs)
  save_rows(out_diffs_full, ratios, timediffs)
  save_columns(out_diffs_avgs,
               ratios,
               mean.(timediffs),
               stderr_of.(timediffs))
end

if isempty(warmup_misses)
  printstyled("Warm-up clean: no compilation inside any timed run.\n"; color = :green)
else
  printstyled("WARNING: compilation landed inside timed runs: \
$(join(("$m: $n" for (m, n) in warmup_misses), ", "))\n"; color = :yellow)
end
