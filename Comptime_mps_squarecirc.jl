using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps
using QuantumClifford

using Statistics
using ProgressMeter
using Strided
using LinearAlgebra

Strided.disable_threads()
nthr = Threads.nthreads()

BLAS.set_num_threads(1)
ITensors.Strided.set_num_threads(1)

function circuit_mps(gates, phases, P; cutoff = 0)
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
    ψ = apply(ψ, g; cutoff = cutoff)
  end

  P_cmps = cmps.PauliSum(P)
  return ITensorMPS.expect(ψ, P_cmps)
end

Nsamples = 100
Ndiv = 7
Nmin = 2
Nmax = 16

cutoff_mps = 1e-12

Nrange = Int.(round.(logrange(Nmin, Nmax, length = Ndiv)))
obs_string = "Z₁"

out_full = "output/comptimes_squarecirc_mps_full.txt"
out_avgs = "output/comptimes_squarecirc_mps_avgs.txt"
param_info = Dict(
  "Nrange" => Nrange,
  "Nsamples" => Nsamples,
  "cutoff_mps" => cutoff_mps,
  "obs" => obs_string)
obsname = "⟨$obs_string⟩ computation time (s)"
initialize_output(out_full, "$obsname, $Nsamples samples", param_info)
initialize_output(out_avgs, "$obsname averages", param_info)

printstyled("Getting MPS computation times of ⟨$obs_string⟩ for square circuits, \
N = $Nrange with cutoff $cutoff_mps.\n\
Nsamples = $Nsamples, $nthr threads.\n"; color = :cyan)
prog = Progress(Ndiv * Nsamples; desc = "Computing…")

times_mps = []

for N in Nrange
  P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
  times_N_mps = Float64[]
  for i in 1:Nsamples
    gates, phases = rotation_circuit(N, N)
    ev_mps, evtime_mps, _ = @timed circuit_mps(gates, phases, P; cutoff = cutoff_mps)
    push!(times_N_mps, evtime_mps)
    next!(prog, showvalues = [("N",N), ("sample",i)])
  end
  push!(times_mps, times_N_mps)
end

save_rows(out_full, Nrange, times_mps)

times_avg_mps = mean.(times_mps)
times_err_mps = @. std(times_mps)/sqrt(Nsamples)

save_columns(out_avgs,
             Nrange,
             times_avg_mps,
             times_err_mps)
