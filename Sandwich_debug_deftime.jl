using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random

# Random.seed!(42)

input_str = length(ARGS) >= 1 ? ARGS[1] : "XXXXXXXXXXXX"

observable = stringtopauli(input_str)
obs_string = pp.inttostring(observable.term, length(input_str))

N = length(obs_string)
t = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : Int(floor(2.5*N))
ss = collect(0:t)

paulistrings = [pp.inttosymbol(rand(0:4^N-1), N) for _ in 1:t]
angles = 4π*rand(Float64, (t,)) # Rotation angles, i.e. 2* exponential phases

# === TESTING ===
# paulistrings = [[:Y, :Y, :Y],[:X, :Y, :Z]]
# paulistrings = [[:Y, :X, :X],[:X, :Y, :Z]]
# paulistrings = [[:X, :Y],[:X, :Z]]
# angles = [π/4, π/4] # Rotation angles, i.e. 2* exponential phases
# == ===

qinds = [collect(1:N) for _ in paulistrings]
rotations = pp.PauliRotation.(paulistrings, qinds)
@assert t == length(rotations)

for s in ss
  printstyled("DEBUGGING expectation value of $(N)-qubit string $obs_string through $t rotation layers with CAMPS until t = $(s).\n"; color = :magenta)
  
  ψ = cmps.CAMPS(N)
  ψ_evo, k = DisentangleCAMPS.evolve(ψ, s, rotations, angles; showprogress = false)

  leftover_rotations = rotations[s+1:end]
  leftover_phases = angles[s+1:end]


  println("Propagating Paulis over $(length(leftover_rotations)) gates")

  sandwichstrings = pp.propagate(leftover_rotations, observable, leftover_phases)

  println("Converting $(length(sandwichstrings))-term sum…")
  sandwichsum, conversiontime, _... = @timed cmps.PauliSum(sandwichstrings)
  println("Done in $conversiontime s.")
  # === Show sums before & after conversion ===
  # println("===== FROM =====\n", sandwichstrings)
  # println("===== TO =====\n", sandwichsum)
  println("Computing expectation value of $(length(sandwichstrings))-term sum on CAMPS with bond dims $(ψ_evo.mps)…")
  ev, evtime, _... = @timed cmps.expectation(ψ_evo, sandwichsum)
  println("Done in $evtime s.")
  println("⟨$(obs_string)⟩ = $ev")
  # println("Full MPS:\n", ψ_evo.mps.data)
  println("\n")
end