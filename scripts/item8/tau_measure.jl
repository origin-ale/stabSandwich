# Item 8, decision 1: measure the singular values that the `:chi3` layer-end search sees, to
# place the numerical-zero tolerance τ. Run from the repository root:
#
#   julia -t 1 --project=. scripts/item8/tau_measure.jl N p1,p2,... samples out.jls
#
# Replays the `:camps` layer loop of `Comptime_prob_tm.jl` (`campssrc_circuit_dynamics`) with
# the script's defaults: same initial states (seed rng_seed + s), same doping masks
# (seed hash((rng_seed + s, p))), π/4 XXZ doped with 3π/16 on XX-YY, `thl = 1e-12`,
# `:chi3`/`:snake`. At every layer end, before the layer disentangler, it SVDs the identity
# and each of the 19 candidate gates on every bond of dimension ≠ 1 (a copy orthogonalized at
# the bond, untruncated `apply_and_svd`, as the search does), and stores σ/σ_max. These are the
# spectra of the first search on each bond; once a gate is applied, later bonds of the same
# sweep see a changed state, which this script does not record.
#
# Output: a `Serialization` file with `recs` (one per candidate spectrum: p, sample s, layer,
# truncated flag, bond, cand (0 = identity, 1–19 = the gate set), stored bond dim, rel σ),
# `summary` (one per run: doped triples, layer at which `truncated` flipped (0 = never), gates
# applied by the layer disentangler, final max bond dim, seconds) and N. Analyse with
# `tau_analyse.jl`.
#
# Uses CampsPP internals (`DisentangleCAMPS._ID_2Q`, `BondSearchDisentangler(...).matrices`,
# `layer_disentangler`, `layer_disentangle!`); adapt if item 8 renames them.
using CampsPP, CampsPP.DisentangleCAMPS
using ITensors, ITensorMPS, QuantumClifford
import MagicTensors as mt
using LinearAlgebra: diag
using Random: MersenneTwister
using Printf, Serialization

const D = CampsPP.DisentangleCAMPS
length(ARGS) == 4 || error("usage: tau_measure.jl N p1,p2,... samples out.jls")
N = parse(Int, ARGS[1])
probs = parse.(Float64, split(ARGS[2], ","))
nsamples = parse(Int, ARGS[3])
outfile = ARGS[4]
# Comptime_prob_tm.jl defaults
ϕ = θ = π / 4; μ = 0.6; rng_seed = 100; thl = 1e-12; magic_phase = 3π / 16

layer_ends = layerends(N, N ÷ 2, xxz_circuit)
de = CampsPP.layer_disentangler(computationalcamps(N, Int[]), :chi3, :snake, N)
mats = vcat([D._ID_2Q], D.BondSearchDisentangler(:chi3).matrices) # cand 0 = identity

recs = NamedTuple[]
summary = NamedTuple[]
for p in probs, s in 1:nsamples
  onebitinds = domainwallstate(MersenneTwister(rng_seed + s), N, μ)[2]
  gates, phases = xxz_circuit(ϕ, θ, N / 2, N)
  mask = doping_mask(MersenneTwister(hash((rng_seed + s, p))), phases, p)
  phases = xy_magic(phases, mask; magicphase = magic_phase)
  ψ = computationalcamps(N, onebitinds)
  norm0_sq = norm(ψ)^2
  truncated = false; trunc_layer = 0; layer = 1; acted = 0
  t0 = time()
  for i in eachindex(gates)
    apply!(ψ, PauliOperator(getpauli(gates[i], N)), phases[i], thl)
    if !truncated && norm(ψ)^2 < norm0_sq * (1 - thl) # as campssrc_circuit_dynamics
      truncated = true; trunc_layer = layer
    end
    i == layer_ends[layer] || continue
    mps = mt.get_mps(ψ)
    dims = linkdims(mps)
    for b in 1:N-1
      dims[b] == 1 && continue
      trial = copy(mps); orthogonalize!(trial, b)
      for (c, M) in enumerate(mats)
        _, S, _ = mt.apply_and_svd(trial, M, b)
        σ = Float64[x for x in diag(S)]
        push!(recs, (p = p, s = s, layer = layer, truncated = truncated, bond = b, cand = c - 1,
                     dim = dims[b], rel = σ ./ maximum(σ)))
      end
    end
    truncated || (acted += first(CampsPP.layer_disentangle!(ψ, de)))
    layer += 1
  end
  push!(summary, (p = p, s = s, ndoped = count(mask), trunc_layer = trunc_layer, acted = acted,
                  maxdim = maxlinkdim(mt.get_mps(ψ)), secs = time() - t0))
  @printf("N=%d p=%.3f s=%d doped=%d trunc_layer=%d acted=%d maxdim=%d %.1fs\n", N, p, s,
          count(mask), trunc_layer, acted, maxlinkdim(mt.get_mps(ψ)), time() - t0)
  flush(stdout)
end
serialize(outfile, (recs = recs, summary = summary, N = N))
