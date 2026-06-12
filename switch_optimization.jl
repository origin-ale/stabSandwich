using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps

using Statistics
using ProgressMeter
using Strided
using LinearAlgebra

Strided.disable_threads()
nthr = Threads.nthreads()

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

function circuit_sandwich(rotations, phases, P, c)
  N = P.nqubits
  t = length(rotations)
  ψ_evo, _ = evolve(cmps.CAMPS(N), t-c, rotations, phases)

  leftover_rotations, leftover_angles = leftover_rotgates(t-c, rotations, phases)

  sandwichstrings = pp.propagate(leftover_rotations, P, leftover_angles)

  return cmps.expectation(ψ_evo, sandwichstrings)
end

N = 100
M = 12
maxc = 5
Nsamples = 50

t = N+M
switchpoints = collect((M-maxc):M)
observable = pp.PauliString(N, [:Z], [1])
obs_string = "Z₁"

out_full = "output/deep_switch_optimization_full_$(N).txt"
out_avgs = "output/deep_switch_optimization_avgs_$(N).txt"
param_info = Dict(
  "N" => N,
  "M" => M,
  "t" => t,
  "switchpoints" => switchpoints,
  "Nsamples" => Nsamples,
  "obs" => obs_string)
obsname = "⟨$obs_string⟩ computation time (s)"
initialize_output(out_full, "$obsname, $Nsamples samples", param_info)
initialize_output(out_avgs, "$obsname averages", param_info)

printstyled("Characterizing computation time of ⟨$obs_string⟩ on $(N)×$(t) random \
rotation circuit with CAMPS-PP.\n\
Switch points: $switchpoints. Nsamples = $Nsamples, $nthr threads.\n"; color = :cyan)
prog = Progress(length(switchpoints) * Nsamples; desc = "Computing…")

times = []
rotations_precomp, phases_precomp = rotation_circuit(2, N)
_ = @timed circuit_sandwich(rotations_precomp, phases_precomp, observable, 1)

for c in switchpoints
  times_c = Float64[]
  for i in 1:Nsamples
    rotations, phases = rotation_circuit(t, N)
    ev, tc, _ = @timed circuit_sandwich(rotations, phases, observable, c)
    push!(times_c, tc)
    next!(prog, showvalues = [("c",c), ("sample",i)])
  end
  push!(times, times_c)
end

save_rows(out_full, switchpoints, times)

times_avg = mean.(times)
times_err = @. std(times)/sqrt(Nsamples)

save_columns(out_avgs, switchpoints, times_avg, times_err)
