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

N = 40
Nsamples = 50
ratiostart = .7
ratioend = 1.1
ratiostep = .1
Nratios = Int(round((ratioend-ratiostart)/ratiostep))+1
progressbar = Progress(Nratios*Nsamples; desc = "Computing…")

times_camps = []
times_pp = []

println("Getting computation times for $Nratios aspect ratios of N=$N t-rotation circuits")

for ratio in ratiostart:ratiostep:ratioend
  t = Int(round(ratio*N))
  P = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
  times_N_camps = Float64[ratio]
  times_N_pp = Float64[ratio]
  for i in 1:Nsamples
    gates, phases = rotation_circuit(t, N)
    ev_camps, evtime_camps, _ = @timed circuit_camps(gates, phases, P)
    push!(times_N_camps, evtime_camps)
    ev_pp, evtime_pp, _ = @timed circuit_pp(gates, phases, P)
    push!(times_N_pp, evtime_pp)
    next!(progressbar, showvalues = [("t/N",ratio), ("sample",i)])
  end
  push!(times_camps, times_N_camps)
  push!(times_pp, times_N_pp)
end
save([times_camps, times_pp], "output/comptimes_ratios_$(N)_full.txt")

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
ratios = collect(ratiostart:ratiostep:ratioend)

save_five_columns(ratios,
                  times_avg_camps,
                  times_err_camps,
                  times_avg_pp,
                  times_err_pp,
                  "output/comptimes_ratios_$(N)_avgs.txt")

println("")