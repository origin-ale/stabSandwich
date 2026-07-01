using Revise
using CampsPP
using DisentangleCAMPS
using ITensors, ITensorMPS

import PauliPropagation as pp
import CliffordMPS as cmps
using QuantumClifford

using Random: MersenneTwister
using ProgressMeter
using Printf

using Strided
using LinearAlgebra
Strided.disable_threads()
nthr = Threads.nthreads()
BLAS.set_num_threads(nthr)
ITensors.Strided.set_num_threads(nthr)

# == Parameters ===============================================================
magic_prob = 1
magic_mode = :xy # Dope on XX-YY with 3π/16 or on ZZ with π/3
μ = 0.6

N = 14
seed = 1

ϕ = π/4
θ = π/4

# CAMPS
camps_thl = 1e-7 # CAMPS SVD truncation threshold
camps_crit = :chi3 # CAMPS entanglement criterion (:chi3 or :entangle)
camps_strat = :snake # CAMPS disentangler strategy (:full, :brickwork, :snake)

if magic_mode == :xy
  magic_phase = 3π/16
  magic_doping = xy_magic
elseif magic_mode == :z
  magic_phase = 3π/16
  magic_doping = z_magic
else
  error("Unrecognised magic_mode: $magic_mode")
end

# == Initialization (as in Comptime_tm.jl init_camps) =========================
rng = MersenneTwister(seed)
ψ_ext, onebitinds = domainwallstate(rng, N, μ)
tm = transferredmagnetization(N, onebitinds)

gates, phases = xxz_circuit(ϕ, θ, N/2, N)
phases = magic_doping(phases, magic_prob; magicphase = magic_phase)

layer_ends = layerends(N, N/2, xxz_circuit)

obs = tm
criterion = camps_crit
strategy = camps_strat
thl = camps_thl
showprogress = true

ψ = deepcopy(ψ_ext)
N = length(ψ)
M = length(gates)
i = 0
layer = 0
evs_camps = []
progress = ProgressUnknown(desc = "Evolving CAMPS (SVD)… gate ", enabled = showprogress)
obs_cmps = cmps.PauliSum(obs)
crit = cmps.DisentangleCriterion(criterion)
strat = cmps.DisentangleStrategy(strategy)

push!(evs_camps, real(cmps.expectation(ψ, obs_cmps)))
layer += 1

while i < M
  global i += 1
  gate = PauliOperator(getpauli(gates[i], N))
  phase = phases[i]

  apply!(ψ, gate, phase, thl)

  if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
    push!(evs_camps, real(cmps.expectation(ψ, obs_cmps)))
    diffs, _ = cmps.disentangle!(ψ, strat, N; criterion = crit, min_diff = 1e-6)
    global layer += 1
    @show diffs
  end

  bd = DisentangleCAMPS.bonddim(ψ)
  next!(progress; showvalues = () -> [("Bond dimension", bd)])
end

finish!(progress)

println("CAMPS (SVD) stopped at gate $i (end of circuit)")
println("Transferred magnetization per layer:")
for (l, ev) in enumerate(evs_camps)
  println("  layer $(l-1): $ev")
end
