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

Nsamples = 100
Ndiv = 13
Nmin = 2
Nmax = 128
Nmax_pp_run = 50 # Pauli prop. is only run for N ≤ Nmax_pp_run

thl_pp = 1e-10

Nrange = Int.(round.(logrange(Nmin, Nmax, length = Ndiv)))
Nrange_pp = filter(≤(Nmax_pp_run), Nrange)
obs_string = "Z₁"

out_full = "output/comptimes_squarecirc_full.txt"
out_avgs = "output/comptimes_squarecirc_avgs.txt"
param_info = Dict(
  "Nrange" => Nrange,
  "Nsamples" => Nsamples,
  "Nmax_pp_run" => Nmax_pp_run,
  "thl_pp" => thl_pp,
  "obs" => obs_string)
obsname = "⟨$obs_string⟩ computation time (s)"
initialize_output(out_full, "$obsname, $Nsamples samples", param_info)
initialize_output(out_avgs, "$obsname averages", param_info)

printstyled("Getting CAMPS and Pauli prop. computation times of ⟨$obs_string⟩ for \
square circuits, N = $Nrange.\n\
Nsamples = $Nsamples, $nthr threads.\n"; color = :cyan)
prog = Progress(Ndiv * Nsamples; desc = "Computing…")

times_camps = []
times_pp = []

for N in Nrange
  P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
  times_N_camps = Float64[]
  times_N_pp = Float64[]
  for i in 1:Nsamples
    gates, phases = rotation_circuit(N, N)
    ev_camps, evtime_camps, _ = @timed circuit_camps(gates, phases, P)
    push!(times_N_camps, evtime_camps)
    if N ≤ Nmax_pp_run
      ev_pp, evtime_pp, _ = @timed circuit_pp(gates, phases, P, thl_pp)
      push!(times_N_pp, evtime_pp)
    end
    next!(prog, showvalues = [("N",N), ("sample",i)])
  end
  push!(times_camps, times_N_camps)
  if N ≤ Nmax_pp_run
    push!(times_pp, times_N_pp)
  end
end

save_rows(out_full, Nrange, times_camps; blockend = true)
save_rows(out_full, Nrange_pp, times_pp; blockend = true)

times_avg_camps = mean.(times_camps)
times_err_camps = @. std(times_camps)/sqrt(Nsamples)
times_avg_pp = mean.(times_pp)
times_err_pp = @. std(times_pp)/sqrt(Nsamples)

save_columns(out_avgs,
             Nrange,
             times_avg_camps,
             times_err_camps,
             times_avg_pp,
             times_err_pp)
