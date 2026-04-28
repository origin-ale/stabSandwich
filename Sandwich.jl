import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random

Random.seed!(42)

include("CiffordMPSPauliPropagation.jl")
include("Interface.jl")

input_str = length(ARGS) >= 1 ? ARGS[1] : "XXXXXXXXXXXX"
χ = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 128

observable = stringtopauli(input_str)
obs_string = pp.inttostring(observable.term, length(input_str))

N = length(obs_string)
t = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : Int(2.5*N)
ψ = cmps.CAMPS(N)

# make these non-random for testing
paulistrings = [pp.inttosymbol(rand(0:4^N-1), N) for _ in 1:t]
qinds = [shuffle(1:N) for _ in 1:t]
rotations = pp.PauliRotation.(paulistrings, qinds)
phases = 2π*rand(Float64, (t,))

printstyled("Calculating expectation value of $(N)-qubit string $obs_string through $t random rotation layers.\n"; color = :cyan)

ψ_evo, k, tstop = DisentangleCAMPS.evolve_bonddim(ψ, χ, rotations, phases; showprogress = true)

leftover_rotations = rotations[tstop+1:end]
leftover_phases = phases[tstop+1:end]

if tstop == t
  println("Maximum bond dimension $χ not reached before s=t=$(t).")
else
  println("Maximum bond dimension $χ reached at s=$(tstop).")
  println("Propagating Paulis over $(length(leftover_rotations)) gates…")
end

sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_phases)

println("Converting $(length(sandwichstrings))-term sum…")
sandwichsum, conversiontime, _... = @timed cmps.PauliSum(sandwichstrings)
println("Done in $conversiontime s.")
println("Calculating expectation value of $(length(sandwichstrings))-term sum on CAMPS with bond dims $(ψ_evo.mps)…")
ev, evtime, _... = @timed cmps.expectation(ψ_evo, sandwichsum)
println("Done in $evtime s.")
println("⟨$(obs_string)⟩ = $ev")