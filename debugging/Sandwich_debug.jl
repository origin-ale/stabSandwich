using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random
using Revise

N = 6
input_str = join(rand(["I","X","Y","Z"], N))

observable = stringtopauli(input_str)
obs_string = pp.inttostring(observable.term, length(input_str))

N = length(obs_string)
t = Int(floor(2.5*N))
ss = collect(0:t)

paulistrings = [pp.inttosymbol(rand(0:4^N-1), N) for _ in 1:t]
phases = 2π * rand(Float64, (t,)) # Exponential phases, ie. -1/2 * rotation angles

qinds = [collect(1:N) for _ in paulistrings]
rotations = pp.PauliRotation.(paulistrings, qinds)
@assert t == length(rotations)

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

  println("Converting $(length(sandwichstrings))-term sum…")
  sandwichsum, conversiontime, _... = @timed cmps.PauliSum(sandwichstrings)
  println("Done in $conversiontime s.")
  println("Computing expectation value of $(length(sandwichstrings))-term sum on CAMPS with bond dims $(ψ_evo.mps)…")
  ev, evtime, _... = @timed cmps.expectation(ψ_evo, sandwichsum)
  println("Done in $evtime s.")
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