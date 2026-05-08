using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using ProgressMeter
using Random

Random.seed!(42)

χs = [1,2]
t = 2
N = 2
paulistrings = [[:X, :Y], [:Y, :X]]
# paulistrings = [[:I, :X], [:X, :Z], [:Y, :Z], [:I, :Z], [:X, :I],
#                 [:Y, :Z], [:Y, :I], [:Y, :Y], [:I, :I], [:Z, :Z],
#                 [:Z, :I], [:X, :X], [:Y, :X], [:Z, :I], [:Y, :X],
#                 [:X, :Y], [:X, :Z], [:Z, :Z], [:X, :Z], [:X, :Y]]
angles = fill(-π/2, t) # Rotation angles, i.e. 2* exponential phases
# angles = [-π/2, -π/2] # Rotation angles, i.e. 2* exponential phases
qinds = [collect(1:N) for _ in paulistrings]
rotations = pp.PauliRotation.(paulistrings, qinds)
@assert t == length(rotations)

ppres = []
campsres = []

input_strs = ["II", "IX", "IY", "IZ", "XI", "XX", "XY", "XZ", "YI", "YX", "YY", "YZ", "ZI", "ZX", "ZY", "ZZ"]
for input_str in input_strs
  observable = stringtopauli(input_str)
  obs_string = pp.inttostring(observable.term, length(input_str))

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
    # (χ == 1) && @show sandwichstrings

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
    (χ == 1) && push!(ppres, ev)
    (χ == 2) && push!(campsres, ev)
  end
end
println(isapprox.(ppres, campsres; atol = 1e-10))