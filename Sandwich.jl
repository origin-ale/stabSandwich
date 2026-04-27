import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random

include("./CiffordMPSPauliPropagation.jl")

N = 24
t = 3*N
χ = 128
ψ = cmps.CAMPS(N)
paulistrings = [pp.inttosymbol(rand(0:4^N-1), N) for _ in 1:t]
qinds = [shuffle(1:N) for _ in 1:t]
rotations = pp.PauliRotation.(paulistrings, qinds)
phases = 2π*rand(Float64, (t,))

println("$(N)-qubit circuit with $t layers.")

ψ_evo, k, tstop = DisentangleCAMPS.evolve_bonddim(ψ, χ, rotations, phases; showprogress = true)

leftover_rotations = rotations[tstop+1:end]
leftover_phases = phases[tstop+1:end]
observable = pp.PauliString(N, [:X, :Y, :Z], [2, N÷2, N-1])

if tstop == t
  println("Maximum bond dimension $χ not reached before s=t=$(t).")
else
  println("Maximum bond dimension $χ reached at s=$(tstop).")
  println("Propagating Paulis over $(length(leftover_rotations)) gates…")
end

sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_phases)

print("Need to sample $(length(sandwichstrings)) strings on $ψ")