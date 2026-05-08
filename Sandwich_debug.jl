using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random
using Revise

# Random.seed!(42)

# input_str = "XXXXXXXXXXXX"
# χs = [4, 8, 16, 32, 64, 128]

input_str = "ZZZZZZZ"
χs = [4, 8, 16, 32, 64]

observable = stringtopauli(input_str)
obs_string = pp.inttostring(observable.term, length(input_str))

N = length(obs_string)
t = Int(floor(2.5*N))

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

for χ in χs
  printstyled("DEBUGGING expectation value of $(N)-qubit string $obs_string through $t rotation layers with CAMPS bond dim $χ.\n"; color = :magenta)
  # printstyled("DEBUGGING expectation value of $(N)-qubit string $obs_string through $t random rotation layers with CAMPS bond dim $χ.\n"; color = :magenta)
  
  ψ = cmps.CAMPS(N)
  ψ_evo, k, tstop = DisentangleCAMPS.evolve_bonddim(ψ, χ, rotations, angles; showprogress = true)

  leftover_rotations = rotations[tstop+1:end]
  leftover_phases = angles[tstop+1:end]

  if tstop == t
    println("Maximum bond dimension $χ not reached before s=t=$(t).")
  else
    println("Maximum bond dimension $χ reached at s=$(tstop).")
    println("Propagating Paulis over $(length(leftover_rotations)) gates")
    # === Print leftover rotations ===
    # for i in eachindex(leftover_rotations)
    #    println(leftover_rotations[i].symbols, '\t', leftover_phases[i])
    # end
  end

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