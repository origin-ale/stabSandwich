using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random
using Revise

# Random.seed!(42)

N = 12
observable = random_paulistring(N)
obs_string = pp.inttostring(observable.term, N)

χ = 128

t = Int(floor(2.5N))
rotations, phases = rotation_circuit(t, N)

ψ = cmps.CAMPS(N)

printstyled("Calculating expectation value of $(N)-qubit string $obs_string through $t random rotation layers.\n"; color = :cyan)

ψ_evo, k, tstop = DisentangleCAMPS.evolve_bonddim(ψ, χ, rotations, phases; showprogress = true)

leftover_rotations, leftover_angles = leftover_rotgates(tstop, rotations, phases)

if tstop == t
  println("Maximum bond dimension $χ not reached before s=t=$(t).")
else
  println("Maximum bond dimension $χ reached at s=$(tstop).")
  println("Propagating Paulis over $(length(leftover_rotations)) gates…")
end

sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_angles)

ev = cmps.expectation(ψ, sandwichstrings; verbose = true)
println("⟨$(obs_string)⟩ = $ev")