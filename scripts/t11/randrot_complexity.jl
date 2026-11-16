# T11 (optional): `Randrot_complexity.jl` on the old and the new stack, for the
# side-by-side comparison in `compare_randrot.jl`. The computation is the script's
# (lines 36-72 at 121384e): for t/N = 0.1:0.1:1.1, `n_samples` random states built from
# t layers of (random Clifford, then a rotation about a random Pauli string at angle
# 4π·rand), recording max_b S_E(b)/N of the CAMPS's MPS part and the N²-sample estimate
# of M₂/N. The script's plotting is left out (not needed for the comparison).
#
#   old stack, from the pre-port worktree (its own Manifest):
#     julia --project=<stabSandwich-old> scripts/t11/randrot_complexity.jl stack=old out=<file>
#   new stack, from this repository:
#     julia --project=. scripts/t11/randrot_complexity.jl stack=new out=<file>
#   optional: N=12 samples=100 seed=1100
#
# Old stack: `CAMPS(N)`, `evolve_deepcliffords` (2N²-deep random Clifford circuit via
# `QOp(:randCliffCircuit)`), `eEntropys!`, `sEntropy`, all drawing from the global RNG.
# New stack: `computationalcamps(mt.QubitCAMPS, N, [])`, `evolve_deepcliffords(rng, …)`
# (uniformly random Clifford, `random_clifford`), `mt.entanglement_entropy`,
# `mt.stabilizer_entropy` (their Phase 6 replacements, pinned to the same spec in
# TestEntropies.jl). The only change to the script's computation is seeding: the global
# RNG is seeded with `seed` (the script is unseeded), and on the new stack the Cliffords
# come from Xoshiro(seed + 1).
#
# The MPS entanglement depends on the Clifford frame, which the two stacks handle
# differently, so its curves are compared side by side only. M₂ is frame-invariant.

args = Dict(String(kv[1]) => String(kv[2]) for kv in (split(a, "="; limit = 2) for a in ARGS))
issubset(keys(args), ("stack", "out", "N", "samples", "seed")) || error("Unknown argument in $ARGS")
stack = args["stack"]
stack in ("old", "new") || error("stack must be old or new")
out = args["out"]
N = parse(Int, get(args, "N", "12"))
n_samples = parse(Int, get(args, "samples", "100"))
seed = parse(Int, get(args, "seed", "1100"))
fractions = 0.1:0.1:1.1

using QuantumClifford
using ITensors, ITensorMPS
using Statistics
using Random
using Dates

if stack == "old"
  using DisentangleCAMPS
  using CliffordMPS
else
  using CampsPP
  using CampsPP.DisentangleCAMPS: evolve_deepcliffords
  import MagicTensors as mt
end

logline(msg) = (println("[", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "] ", msg); flush(stdout))

Random.seed!(seed)
rng_cliff = Xoshiro(seed + 1)   # new stack only

""" One sample: (max bond entropy / N, sampled M₂ / N) after t layers. """
function sample_old(N, t)
  ψ = CAMPS(N)
  paulistrings = [PauliOperator(0x0, rand(Bool, N), rand(Bool, N)) for _ in 1:t]
  angles = 4π * rand(Float64, (t,))
  ψ_evo, _ = evolve_deepcliffords(ψ, t, paulistrings, angles)
  return maximum(eEntropys!(ψ_evo.mps)) / N, sEntropy(ψ_evo.mps, N^2; α = 2) / N
end

function sample_new(N, t)
  ψ = computationalcamps(mt.QubitCAMPS, N, Int[])
  paulistrings = [PauliOperator(0x0, rand(Bool, N), rand(Bool, N)) for _ in 1:t]
  angles = 4π * rand(Float64, (t,))
  ψ_evo, _ = evolve_deepcliffords(rng_cliff, ψ, t, paulistrings, angles)
  mps = mt.get_mps(ψ_evo)
  return maximum(mt.entanglement_entropy(mps)) / N, mt.stabilizer_entropy(mps, N^2; α = 2) / N
end

sample_fn = stack == "old" ? sample_old : sample_new
logline("Randrot complexity, stack=$stack, N=$N, samples=$n_samples, seed=$seed, Julia $VERSION")

open(out, "w") do io
  println(io, "# stack=$stack N=$N samples=$n_samples seed=$seed")
  println(io, "# t/N avg_EE err_EE avg_SRE err_SRE time_s")
end
for f in fractions
  t = round(Int, f * N)
  ee, sre = Float64[], Float64[]
  time = @elapsed for _ in 1:n_samples
    e, m = sample_fn(N, t)
    push!(ee, e); push!(sre, m)
  end
  se(x) = std(x) / sqrt(length(x))
  open(out, "a") do io
    println(io, "$f $(mean(ee)) $(se(ee)) $(mean(sre)) $(se(sre)) $time")
  end
  logline("t/N=$f (t=$t): EE/N $(round(mean(ee), digits = 4)), M2/N $(round(mean(sre), digits = 4)), $(round(time, digits = 1)) s")
end
logline("Done; results in $out")
