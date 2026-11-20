# Item 8 Verification: accuracy of `:camps` against dense, at small N. Run from the repository
# root, with --project set to the tree under test (this tree, or a checkout of the commit
# before item 8):
#
#   julia -t 1 --project=<tree> scripts/item8/accuracy_dense.jl N p1,p2,... samples out.jls [thl] [mode]
#
# Same samples as `Comptime_prob_tm.jl` and `tau_measure.jl`: initial state seed rng_seed + s,
# doping mask seed hash((rng_seed + s, p)), π/4 XXZ doped with 3π/16 on XX-YY, N/2 layers,
# `thl` (default 1e-12, as the script), `:chi3`/`:snake`, observable the transferred
# magnetization. Each sample runs the CAMPS dynamics and the dense reference
# (CampsPP/test/DisentangleCAMPS/DenseReference.jl, loaded from the tree this script is in) at
# the layer ends.
#
# mode (default `campssrc`): `campssrc` calls `campssrc_circuit_dynamics` (track = true).
# `noflag` replays its loop (`apply!`, record ⟨O⟩, layer disentangler, bond dimension) but
# always disentangles, ignoring the `truncated` flag. It also records whether that flag would
# have tripped and the number of gates the disentangler applied. Use it with thl = 0, where
# rounding alone can trip the flag and switch the disentangler off (see the comment in
# `campssrc_circuit_dynamics`). It uses CampsPP internals (`convert_observable`,
# `expectation_value`, `layer_disentangler`, `layer_disentangle!`).
#
# Output: a `Serialization` file with `recs` (one per sample: p, s, doped triples, |⟨O⟩ − dense|
# per record point, max bond dimension per layer end, seconds; with `noflag` also `acted` and
# `flag`) and N, thl, mode. Compare two files with `accuracy_compare.jl`.
using CampsPP, CampsPP.DisentangleCAMPS
using QuantumClifford
import PauliPropagation as pp
using ITensors, ITensorMPS
import MagicTensors as mt
using Random: MersenneTwister
using Printf, Serialization

include(joinpath(@__DIR__, "..", "..", "CampsPP", "test", "DisentangleCAMPS", "DenseReference.jl"))
const DR = DenseReference

4 <= length(ARGS) <= 6 || error("usage: accuracy_dense.jl N p1,p2,... samples out.jls [thl] [campssrc|noflag]")
N = parse(Int, ARGS[1])
probs = parse.(Float64, split(ARGS[2], ","))
nsamples = parse(Int, ARGS[3])
outfile = ARGS[4]
thl = length(ARGS) >= 5 ? parse(Float64, ARGS[5]) : 1e-12
mode = length(ARGS) >= 6 ? Symbol(ARGS[6]) : :campssrc
mode in (:campssrc, :noflag) || error("mode must be campssrc or noflag")
N <= 14 || error("dense reference: N ≤ 14")
# Comptime_prob_tm.jl defaults
ϕ = θ = π / 4; μ = 0.6; rng_seed = 100; magic_phase = 3π / 16

# As CampsPP/test/helpers.jl (the dense reference's string conventions, CONVENTIONS.md §2, §3)
function gatestring(rot, N)
  letters = fill('I', N)
  for (sym, q) in zip(rot.symbols, rot.qinds)
    letters[q] = string(sym)[1]
  end
  return String(letters)
end
termstring(term, N) = String(["IXYZ"[Int((term >> (2 * (j - 1))) & 0x3) + 1] for j in 1:N])
dense_expect(ψ, O::pp.PauliSum) =
  sum(coeff * DR.expect(ψ, termstring(term, O.nqubits)) for (term, coeff) in O.terms)
function dense_dynamics(N, ones, gates, phases, O, record_at)
  ψ = DR.basis_state(N, ones)
  evs = [dense_expect(ψ, O)]
  for (i, (g, φ)) in enumerate(zip(gates, phases))
    DR.rotate!(ψ, gatestring(g, N), φ)
    i in record_at && push!(evs, dense_expect(ψ, O))
  end
  return evs
end

# The loop of `campssrc_circuit_dynamics`, always disentangling. Returns (evs, bds, acted, flag).
function camps_noflag(ψ0, gates, phases, thl, O, layer_ends)
  ψ = deepcopy(ψ0)
  N = length(mt.get_mps(ψ))
  obs = CampsPP.convert_observable(ψ, O)
  de = CampsPP.layer_disentangler(ψ, :chi3, :snake, N)
  norm0_sq = norm(ψ)^2
  flag = false; acted = 0; layer = 1
  evs = [real(CampsPP.expectation_value(ψ, obs))]
  bds = [DisentangleCAMPS.bonddim(ψ)] # t = 0, as campssrc_circuit_dynamics with track = true
  for i in eachindex(gates)
    apply!(ψ, PauliOperator(getpauli(gates[i], N)), phases[i], thl)
    flag |= norm(ψ)^2 < norm0_sq * (1 - thl)
    layer <= length(layer_ends) && i == layer_ends[layer] || continue
    push!(evs, real(CampsPP.expectation_value(ψ, obs)))
    acted += first(CampsPP.layer_disentangle!(ψ, de))
    push!(bds, DisentangleCAMPS.bonddim(ψ))
    layer += 1
  end
  return evs, bds, acted, flag
end

layer_ends = layerends(N, N ÷ 2, xxz_circuit)
recs = NamedTuple[]
for p in probs, s in 1:nsamples
  onebitinds = domainwallstate(MersenneTwister(rng_seed + s), N, μ)[2]
  gates, phases = xxz_circuit(ϕ, θ, N / 2, N)
  mask = doping_mask(MersenneTwister(hash((rng_seed + s, p))), phases, p)
  phases = xy_magic(phases, mask; magicphase = magic_phase)
  O = transferredmagnetization(N, onebitinds)
  t0 = time()
  ψ0 = computationalcamps(N, onebitinds)
  acted = flag = missing
  if mode == :campssrc
    _, _, evs, bds = campssrc_circuit_dynamics(ψ0, gates, phases, thl, O; criterion = :chi3,
      strategy = :snake, layer_ends = layer_ends, track = true)
  else
    evs, bds, acted, flag = camps_noflag(ψ0, gates, phases, thl, O, layer_ends)
  end
  secs = time() - t0
  ref = dense_dynamics(N, onebitinds, gates, phases, O, layer_ends)
  length(ref) == length(evs) || error("record points differ: $(length(ref)) vs $(length(evs))")
  err = abs.(Float64.(real.(evs)) .- real.(ref))
  push!(recs, (p = p, s = s, ndoped = count(mask), err = err, bds = bds, secs = secs,
               acted = acted, flag = flag))
  @printf("N=%d p=%.3f s=%d doped=%d max|Δ|=%.2e maxbd=%d acted=%s flag=%s %.1fs\n", N, p, s,
          count(mask), maximum(err), maximum(bds), acted, flag, secs)
  flush(stdout)
end
serialize(outfile, (recs = recs, N = N, thl = thl, mode = mode))
