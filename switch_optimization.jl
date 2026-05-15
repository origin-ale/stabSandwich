using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random
using Revise

using ProgressMeter
using Statistics

# Random.seed!(42)

function circuit_sandwich(rotations, phases, P, c)
  N = P.nqubits
  t = length(rotations)
  ψ_evo, k = evolve(cmps.CAMPS(N), t-c, rotations, phases)

  leftover_rotations, leftover_angles = leftover_rotgates(t-c, rotations, phases)

  sandwichstrings = pp.propagate(leftover_rotations, P, leftover_angles)

  return cmps.expectation(ψ_evo, sandwichstrings)
end

function save_rows(filename, params, vals)
  open(filename, "w") do f
    for i in eachindex(vals)
      print(f,"$(params[i]) ")
      for e in vals[i]
        print(f, "$e ")
      end
      println(f, "")
    end
  end
end

function save_columns(filename, params, vals...)
  open(filename, "w") do f
    for i in eachindex(params)
      print(f,"$(params[i]) ")
      for e in vals
        print(f, "$(e[i]) ")
      end
      println(f, "")
    end
  end
end

N = 60
switchpoints = 0:1:12
Nsamples = 50
out_full = "output/switch_optimization_full_$(N).txt"
out_avgs = "output/switch_optimization_avgs_$(N).txt"

t = N
observable = pp.PauliString(N, [:Z], [1])
obs_string = "Z₁"
prog = Progress(length(switchpoints)*Nsamples; desc = "Computing…" )

printstyled("Characterizing computation time of ⟨$(obs_string)⟩ on $(N)×$(t) random rotation circuit with CAMPS-PP.
Switch points: ", collect(switchpoints), "\n"; color = :cyan)

times = []
rotations_precomp, phases_precomp = rotation_circuit(2, N)
_ = @timed circuit_sandwich(rotations_precomp, phases_precomp, observable, 1)

for c in switchpoints
  times_c = []
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