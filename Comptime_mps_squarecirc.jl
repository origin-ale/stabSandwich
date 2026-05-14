using Statistics

using ITensors, ITensorMPS
using CampsPP
using DisentangleCAMPS
using ProgressMeter

import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford

function save(arr, filename)
  open(filename, "w") do f
    for row in arr
      println(f, join(row, " "))
    end
  end
end

function save_three_columns(a, b, c, filename)
  n = max(length(a), length(b), length(c))
  open(filename, "w") do f
    for i in 1:n
      ai = i <= length(a) ? a[i] : ""
      bi = i <= length(b) ? b[i] : ""
      ci = i <= length(c) ? c[i] : ""
      println(f, "$(ai) $(bi) $(ci)")
    end
  end
end

function circuit_mps(gates, phases, P; cutoff = 0)
  N = P.nqubits
  I = PauliOperator(0x0, fill(false,N), fill(false,N))

  gate_paulis = cmps.PauliOperator.(getpauli.(gates, N))

  sites = siteinds("Qubit", N)
  states = ["Up" for i in sites]
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
cutoff = 1e-12
progressbar = Progress(Ndiv*Nsamples; desc = "Computing…")

times_mps = []
Nrange = Int.(round.(logrange(2, 16, length=Ndiv)))

println("Getting MPS computation times for square circuits, N = $Nrange with cutoff $cutoff")

for N in Nrange
  P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
  times_N_mps = Float64[N]
  for i in 1:Nsamples 
    gates, phases = rotation_circuit(N, N)
    ev_mps, evtime_mps, _ = @timed circuit_mps(gates, phases, P; cutoff = cutoff)
    push!(times_N_mps, evtime_mps)
    next!(progressbar, showvalues = [("N",N), ("sample",i)])
  end
  push!(times_mps, times_N_mps)
end
save(times_mps, "output/comptimes_squarecirc_mps_full.txt")

for arr in times_mps
  popat!(arr, 1)
end

times_avg_mps = mean.(times_mps)
times_err_mps = @. std(times_mps)/sqrt(Nsamples)

save_three_columns(Nrange,
                  times_avg_mps,
                  times_err_mps,
                  "output/comptimes_squarecirc_mps_avgs.txt")