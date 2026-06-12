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

function circuit_camps(gates, phases, P)
  ψ, _ = evolve(cmps.CAMPS(P.nqubits), length(phases), gates, phases)
  return cmps.expectation(ψ, P)
end

function circuit_pp(gates, phases, P, thl)
  propagated = pp.propagate(gates, P, -2 .* phases; min_abs_coeff = thl)
  return pp.overlapwithzero(propagated)
end

N = 40
Nsamples = 50
ratiostart = .7
ratioend = 1.1
ratiostep = .1

thl_pp = 1e-10

ratios = collect(ratiostart:ratiostep:ratioend)
Nratios = length(ratios)
P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
obs_string = "Z₁"

out_full = "output/comptimes_ratios_$(N)_full.txt"
out_avgs = "output/comptimes_ratios_$(N)_avgs.txt"
param_info = Dict(
  "N" => N,
  "Nsamples" => Nsamples,
  "ratios" => ratios,
  "thl_pp" => thl_pp,
  "obs" => obs_string)
obsname = "⟨$obs_string⟩ computation time (s)"
initialize_output(out_full, "$obsname, $Nsamples samples", param_info)
initialize_output(out_avgs, "$obsname averages", param_info)

printstyled("Getting CAMPS and Pauli prop. computation times of ⟨$obs_string⟩ for \
$Nratios aspect ratios of N=$N t-rotation circuits.\n\
Nsamples = $Nsamples, $nthr threads.\n"; color = :cyan)
prog = Progress(Nratios * Nsamples; desc = "Computing…")

times_camps = []
times_pp = []

for ratio in ratios
  t = Int(round(ratio*N))
  times_t_camps = Float64[]
  times_t_pp = Float64[]
  for i in 1:Nsamples
    gates, phases = rotation_circuit(t, N)
    ev_camps, evtime_camps, _ = @timed circuit_camps(gates, phases, P)
    push!(times_t_camps, evtime_camps)
    ev_pp, evtime_pp, _ = @timed circuit_pp(gates, phases, P, thl_pp)
    push!(times_t_pp, evtime_pp)
    next!(prog, showvalues = [("t/N",ratio), ("sample",i)])
  end
  push!(times_camps, times_t_camps)
  push!(times_pp, times_t_pp)
end

save_rows(out_full, ratios, times_camps; blockend = true)
save_rows(out_full, ratios, times_pp; blockend = true)

times_avg_camps = mean.(times_camps)
times_err_camps = @. std(times_camps)/sqrt(Nsamples)
times_avg_pp = mean.(times_pp)
times_err_pp = @. std(times_pp)/sqrt(Nsamples)

save_columns(out_avgs,
             ratios,
             times_avg_camps,
             times_err_camps,
             times_avg_pp,
             times_err_pp)
