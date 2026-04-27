import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random

include("./CiffordMPSPauliPropagation.jl")

N = 12
t = N
ψ = cmps.CAMPS(N)
paulistrings = [pp.inttosymbol(rand(0:4^N-1), N) for _ in 1:t]
qinds = [shuffle(1:N) for _ in 1:t]
rotations = pp.PauliRotation.(paulistrings, qinds)
phases = 2π*rand(Float64, (t,))

ψ_evo, k = evolve(ψ, t÷2, rotations, phases; showprogress = true)
leftover_rotations = rotations[t÷2+1:end]
leftover_phases = phases[t÷2+1:end]

observable = pp.PauliString(N, [:X, :Y, :Z], [2, N÷2, N-1])
println("Propagating Paulis over $(length(leftover_rotations)) gates…")
sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_phases)
print("Need to sample $(length(sandwichstrings)) strings on $ψ")