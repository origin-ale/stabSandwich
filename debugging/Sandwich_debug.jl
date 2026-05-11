using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random
using Revise

N = 6
observable = random_paulistring(N)

t = Int(floor(2.5*N))
ss = collect(0:t)
rotations, phases = rotation_circuit(t, N)

printstyled(repeat("=",80), "\n\n"; color = :cyan)
printstyled("DEBUGGING expectation value of $(N)-qubit string $obs_string through $t \
  random rotation layers.\n"; color = :cyan)

evs = []

for s in ss
  printstyled(repeat("-",16), " CAMPS until t = $s ", repeat("-",16), "\n"; color = :light_cyan)
  
  ψ = cmps.CAMPS(N)
  ψ_evo, k = DisentangleCAMPS.evolve(ψ, s, rotations, phases; showprogress = true)

  leftover_rotations, leftover_angles = leftover_rotgates(s, rotations, phases)
  sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_angles)

  ev = cmps.expectation(ψ, sandwichstrings; verbose = false)
  println("⟨$(obs_string)⟩ = $ev")
  println("\n")
  push!(evs, ev)
end

if all(isapprox.(evs[1], evs; atol=1e-10))
  printstyled("All expectation values match.\n"; color = :green)
else
  printstyled("Expectation values DO NOT match:\n"; color = :red)
  print(evs)
end