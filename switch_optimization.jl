using CampsPP
using CampsPP.DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import MagicTensors as mt

using Statistics

@isdefined(ScriptSetup) || include(joinpath(@__DIR__, "scripts", "ScriptSetup.jl"))
using .ScriptSetup

nthr = Threads.nthreads()
setup_threads(1)

function circuit_sandwich(rotations, phases, P, c)
  N = P.nqubits
  t = length(rotations)
  ψ_evo, _ = evolve(computationalcamps(N, Int[]), t-c, rotations, phases)

  leftover_rotations, leftover_angles = leftover_rotgates(t-c, rotations, phases)

  sandwichstrings = pp.propagate(leftover_rotations, P, leftover_angles)

  return mt.expectation(ψ_evo, sandwichstrings)
end

# == Parameters ===============================================================
# Override from the command line as key=value (`help` lists them), e.g.
#   julia --project=. switch_optimization.jl N=20 M=12 Nsamples=50
params = @parameters begin
  N = 20
  M = 12          # circuit depth t = N + M
  maxc = 12       # switch points c = M - maxc, …, M (gates left for PP)
  # Explicit switch points, comma-separated (e.g. switchpoints=0,1,2 for the old square
  # runs, M = 0); replaces the range above when given, so c > M is possible
  switchpoints = Int[]
  Nsamples = 50
  outdir = "output/switch_optimization"
end

t = N+M
isempty(switchpoints) && (switchpoints = collect((M-maxc):M))
all(c -> 0 <= c <= t, switchpoints) ||
  error("Switch points must lie in 0:$t (gates left for PP, out of t = N + M), got $switchpoints")
observable = pp.PauliString(N, [:Z], [1])
obs_string = "Z₁"

mkpath(outdir)
out_full = joinpath(outdir, "M$(M)_switch_optimization_full_$(N).txt")
out_avgs = joinpath(outdir, "M$(M)_switch_optimization_avgs_$(N).txt")
param_info = [param_header(params; exclude = (:outdir, :switchpoints));
  "t" => string(t);
  "switchpoints" => join(switchpoints, ",");
  "obs" => obs_string]
obsname = "⟨$obs_string⟩ computation time (s)"
initialize_output(out_full, "$obsname, $Nsamples samples", param_info)
initialize_output(out_avgs, "$obsname averages", param_info)

log_message("Characterizing computation time of ⟨$obs_string⟩ on $(N)×$(t) random \
rotation circuit with CAMPS-PP. Switch points: $switchpoints. Nsamples = $Nsamples, \
$nthr threads.")
prog = ProgressLog("sample", length(switchpoints) * Nsamples; frac = 0)

times = []
rotations_precomp, phases_precomp = rotation_circuit(2, N)
_ = @timed circuit_sandwich(rotations_precomp, phases_precomp, observable, 1)

for (ic, c) in enumerate(switchpoints)
  times_c = Float64[]
  for i in 1:Nsamples
    rotations, phases = rotation_circuit(t, N)
    ev, tc, _ = @timed circuit_sandwich(rotations, phases, observable, c)
    push!(times_c, tc)
    log_progress!(prog, (ic - 1) * Nsamples + i, "c" => c, "sample" => "$i/$Nsamples")
  end
  push!(times, times_c)

  initialize_output(out_full, "$obsname, $Nsamples samples", param_info)
  save_rows(out_full, switchpoints[1:length(times)], times)

  times_avg = mean.(times)
  times_err = @. std(times)/sqrt(Nsamples)

  initialize_output(out_avgs, "$obsname averages", param_info)
  save_columns(out_avgs, switchpoints, times_avg, times_err)
end
