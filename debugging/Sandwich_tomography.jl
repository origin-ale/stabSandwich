using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Revise

χs = [1,2]
N = 2

paulistrings = [[:X, :Y], [:Y, :X]]
t = length(paulistrings)
phases = fill(π/2, t)
rotations, phases = rotation_circuit(paulistrings, phases, N)

ppres = []
campsres = []

input_strs = ["II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX", "ZY", "ZZ"]
for input_str in input_strs
  observable = stringtopauli_sym(input_str)
  obs_string = pp.inttostring(observable.term, length(input_str))

  for χ in χs
    printstyled("DEBUGGING expectation value of $(N)-qubit string $obs_string through $t rotation layers with CAMPS bond dim $χ.\n"; color = :magenta)
    
    ψ = cmps.CAMPS(N)
    ψ_evo, k, tstop = DisentangleCAMPS.evolve_bonddim(ψ, χ, rotations, phases; showprogress = true)

    leftover_rotations, leftover_angles = leftover_rotgates(tstop, rotations, phases)

    if tstop == t
      println("Maximum bond dimension $χ not reached before s=t=$(t).")
    else
      println("Maximum bond dimension $χ reached at s=$(tstop).")
      println("Propagating Paulis over $(length(leftover_rotations)) gates")
    end

    sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_angles)

    ev = cmps.expectation(ψ, sandwichstrings; verbose = true)
    println("⟨$(obs_string)⟩ = $ev")
    println("\n")
    (χ == 1) && push!(ppres, ev)
    (χ == 2) && push!(campsres, ev)
  end
end

println("Coincidence of traces:")
println(isapprox.(ppres, campsres; atol = 1e-10))