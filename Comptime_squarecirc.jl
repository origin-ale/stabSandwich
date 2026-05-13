# using Revise # Use at own risk, plays bad with world age
using ProgressMeter
using Statistics

using DisentangleCAMPS
using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp

function save(arrays, filename)
  open(filename, "w") do f
    for arr in arrays
      for row in arr
        println(f, join(row, " "))
      end
      println(f)
    end
  end
end

function save_five_columns(a, b, c, d, e, filename)
  n = max(length(a), length(b), length(c), length(d), length(e))
  open(filename, "w") do f
    for i in 1:n
      ai = i <= length(a) ? a[i] : ""
      bi = i <= length(b) ? b[i] : ""
      ci = i <= length(c) ? c[i] : ""
      di = i <= length(d) ? d[i] : ""
      ei = i <= length(e) ? e[i] : ""
      println(f, "$(ai) $(bi) $(ci) $(di) $(ei)")
    end
  end
end

function circuit_camps(gates, phases, P)
  ψ, k = evolve(cmps.CAMPS(P.nqubits), length(phases), gates, phases)
  return cmps.expectation(ψ, P)
end

function circuit_pp(gates, phases, P)
  propagated = pp.propagate(gates, P, -2 .* phases; min_abs_coeff = 1e-10)
  return pp.overlapwithzero(propagated)
end

Nsamples = 100
Ndiv = 13
progressbar = Progress(Ndiv*Nsamples; desc = "Computing…")

times_camps = []
times_pp = []

for N in Int.(ceil.(logrange(2, 128, length=Ndiv)))
  P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
  times_N_camps = Float64[N]
  times_N_pp = Float64[N]
  for i in 1:Nsamples
    gates, phases = rotation_circuit(N, N)
    ev_camps, evtime_camps, _ = @timed circuit_camps(gates, phases, P)
    push!(times_N_camps, evtime_camps)
    if N ≤ 50 # run like this
      ev_pp, evtime_pp, _ = @timed circuit_pp(gates, phases, P)
      push!(times_N_pp, evtime_pp)
    end
    next!(progressbar, showvalues = [("N",N), ("sample",i)])
  end
  push!(times_camps, times_N_camps)
  push!(times_pp, times_N_pp)
end
save([times_camps, times_pp], "output/comptimes_squarecirc_full.txt")

for arr in times_camps
  popat!(arr, 1)
end
for arr in times_pp
  popat!(arr, 1)
end

times_avg_camps = mean.(times_camps)
times_err_camps = @. std(times_camps)/sqrt(Nsamples)
times_avg_pp = mean.(times_pp)
times_err_pp = @. std(times_pp)/sqrt(Nsamples)
sizes = Int.(ceil.(logrange(2, 128, length=Ndiv)))

save_five_columns(sizes,
                  times_avg_camps,
                  times_err_camps,
                  times_avg_pp[1:10],
                  times_err_pp[1:10],
                  "output/comptimes_squarecirc_avgs.txt")