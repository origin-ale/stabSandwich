# Phase 5: local micro-benchmarks of the NEW stack (`main`): the pp → `QubitPauliSum`
# converter, the frame conjugation, the `expectation` seam (conjugate → OpSum → MPO → inner) and
# a per-term prototype of ⟨ψ|O|ψ⟩ (the plan's optional replacement for the MPO route; benchmark
# only, not in src/). Run from the repository root, under a timeout:
#
#   julia -t 1 --project=. scripts/phase5/bench_new.jl [maxmpo=1161] [budget=900] [out=output/phase5/bench_new.tsv]
#
#   maxmpo  largest sum evaluated through the MPO route (it grows roughly as terms^2.75)
#   budget  seconds allowed per per-term run; a run whose estimate (from its first terms)
#           exceeds it is skipped and only the estimate is reported
#
# State: the T8 hand-off state (p = 0.06 sample 1 at its switch gate, bond dimension 112),
# rebuilt as in TestMagicTensorsPP.jl. Sums: the T8 sums plus output/phase5/sums/ (gen_sums.jl).
# One thread, BLAS on one thread. BenchmarkTools comes from scripts/phase5/benchenv.

push!(LOAD_PATH, joinpath(@__DIR__, "benchenv"))
using BenchmarkTools
using CampsPP
import MagicTensors as mt
using ITensors, ITensorMPS
using QuantumClifford: PauliOperator
const DisentangleCAMPS = CampsPP.DisentangleCAMPS
using LinearAlgebra, Strided
Strided.disable_threads()
BLAS.set_num_threads(1)
include("common.jl")

args = Dict(split(a, "="; limit = 2) |> kv -> kv[1] => kv[2] for a in ARGS)
for k in keys(args)
  k in ("maxmpo", "budget", "out") || error("Unrecognised argument `$k` (expected maxmpo, budget, out)")
end
maxmpo = parse(Int, get(args, "maxmpo", "1161"))
budget = parse(Float64, get(args, "budget", "900"))
outpath = get(args, "out", joinpath(@__DIR__, "..", "..", "output", "phase5", "bench_new.tsv"))
flushprintln("threads: julia $(Threads.nthreads()), BLAS $(BLAS.get_num_threads()); maxmpo = $maxmpo, budget = $budget s")

# == State ======================================================================
st = read_state()
N = st.N
gates, phases = sample_circuit(st, xxz_circuit)
t_state = @elapsed (ψ, k, s_rebuilt) = DisentangleCAMPS.evolve_bonddim(
  computationalcamps(mt.QubitCAMPS, N, st.excited), CHI, gates, phases)
s_rebuilt == st.s || error("rebuild switched at gate $s_rebuilt, expected $(st.s)")
bd = DisentangleCAMPS.bonddim(ψ)
flushprintln("state: s = $s_rebuilt, k = $k, bond dimension $bd ($(round(t_state; digits = 1)) s)")

# == Per-term prototype =========================================================
# ⟨ψ|O|ψ⟩ = Σ_P c_P ⟨mps|P|mps⟩ with O conjugated into the MPS frame (the same `mt.conjugate`
# and the same xz → I/X/Y/Z rule as MagicTensors' OpSum route). The MPS is put in canonical
# form with centre at site 1, so the sites to the right of a term's last non-identity site
# contract to the identity and each term costs only up to its last non-identity site.
const OPNAMES = ("Id", "X", "Y", "Z")
xzcode(x::Bool, z::Bool) = x ? (z ? 2 : 1) : (z ? 3 : 0)   # MagicTensors' _int_from_xz

struct PerTermCtx
  ket::Vector{ITensor}
  bra::Vector{ITensor}
  ops::Vector{NTuple{4, ITensor}}
  links::Vector{Index}
  normsq::ComplexF64
end

function PerTermCtx(ψ)
  mps = orthogonalize(mt.get_mps(ψ), 1)
  sites = siteinds(mps)
  ket = [mps[j] for j in 1:N]
  bra = [dag(prime(mps[j])) for j in 1:N]
  ops = [Tuple(op(n, sites[j]) for n in OPNAMES) for j in 1:N]
  links = [linkind(mps, j) for j in 1:N-1]
  return PerTermCtx(ket, bra, ops, links, scalar(ket[1] * dag(ket[1])))
end

"Conjugate into the frame and list (codes, coeff) per term; codes[j] ∈ 0:3 (I, X, Y, Z)."
function frame_terms(ψ, Q::mt.QubitPauliSum)
  Qc = mt.conjugate(Q, mt.get_clifford(ψ; daggered = true))
  out = Vector{Tuple{Vector{Int8}, ComplexF64}}(undef, length(Qc.paulis))
  for (i, (xz, c)) in enumerate(Qc.paulis)
    P = PauliOperator(0x00, N, xz)
    out[i] = ([Int8(xzcode(P[j]...)) for j in 1:N], c)
  end
  return out
end

lastsite(codes) = something(findlast(!iszero, codes), 0)

function close_env(ctx, L, m)
  m == 0 && return ctx.normsq
  m == N && return scalar(L)
  l = ctx.links[m]
  return scalar(L * delta(l, prime(l)))
end

"Naive: every term contracted from site 1 to its last non-identity site."
function perterm_naive(ctx, terms)
  ev = zero(ComplexF64)
  for (codes, c) in terms
    m = lastsite(codes)
    L = ITensor(1.0)
    for j in 1:m
      L = ((L * ctx.ket[j]) * ctx.ops[j][codes[j] + 1]) * ctx.bra[j]
    end
    ev += c * close_env(ctx, L, m)
  end
  return ev
end

"Prefix-shared: terms sorted by their code strings, left environments reused along common prefixes."
function perterm_prefix(ctx, terms)
  order = sortperm(terms; by = first)
  envs = Vector{ITensor}(undef, N)
  prev = fill(Int8(-1), N)
  valid = 0                                    # envs[1:valid] match prev[1:valid]
  ev = zero(ComplexF64)
  for i in order
    codes, c = terms[i]
    m = lastsite(codes)
    j0 = 0
    while j0 < min(valid, m) && codes[j0 + 1] == prev[j0 + 1]
      j0 += 1
    end
    L = j0 == 0 ? ITensor(1.0) : envs[j0]
    for j in j0+1:m
      L = ((L * ctx.ket[j]) * ctx.ops[j][codes[j] + 1]) * ctx.bra[j]
      envs[j] = L
      prev[j] = codes[j]
    end
    valid = m                                  # entries past m belong to an older prefix
    ev += c * close_env(ctx, L, m)
  end
  return ev
end

ctx = PerTermCtx(ψ)

# == Benchmarks =================================================================
med(b) = median(b).time / 1e9
fmt(x) = x isa Number && isfinite(x) ? @sprintf("%.4g", x) : "—"

"Time f on the first `n0` terms, extrapolate to all; run the full benchmark only within budget."
function timed_perterm(f, terms; n0 = 200)
  n = length(terms)
  sub = terms[1:min(n0, n)]
  f(ctx, sub)                                  # compile
  t0 = @elapsed f(ctx, sub)
  est = t0 * n / length(sub)
  est > budget && return (NaN, NaN + 0im, est)
  if est > 60                                  # long: one timed run
    t = @elapsed v = f(ctx, terms)
    return (t, v, est)
  end
  b = @benchmark $f($ctx, $terms) samples = 5 evals = 1 seconds = max(1.0, min(budget, 3 * est))
  return (med(b), f(ctx, terms), est)
end

files = sum_files()
# Warm-up on the smallest sum, so no row below carries compilation.
let S = read_sum(first(files))[1]
  mt.expectation(ψ, S); frame_terms(ψ, mt.QubitPauliSum(S))
end

cols = ["nterms", "conv_s", "conj_s", "seam_s", "opsum_s", "mpo_build_s", "inner_s", "mpo_maxdim",
        "frame_terms_s", "naive_s", "naive_est_s", "prefix_s", "prefix_est_s",
        "ev_seam", "ev_naive", "ev_prefix", "ev_old"]
rows = []
for f in files
  S, header, _ = read_sum(f)
  n = length(S)
  flushprintln("\n== $(basename(f)): $n terms")
  ev_old = haskey(header, "old_expectation_mpo") ?
           parse(ComplexF64, replace(header["old_expectation_mpo"], " " => "")) : NaN + 0im

  bconv = @benchmark mt.QubitPauliSum($S) seconds = 10
  Q = mt.QubitPauliSum(S)
  C = mt.get_clifford(ψ; daggered = true)
  bconj = @benchmark mt.conjugate($Q, $C) seconds = 10
  flushprintln("converter $(fmt(med(bconv))) s, conjugation $(fmt(med(bconj))) s")

  seam = opsum = build = inn = NaN; maxdim_mpo = 0; ev_seam = NaN + 0im
  if n <= maxmpo
    # The seam end to end (one timed run: seconds to minutes, compiled by the warm-up), then
    # one run split into its stages.
    GC.gc()
    seam = @elapsed ev_seam = mt.expectation(ψ, S)
    Qc = mt.conjugate(Q, C)
    opsum = @elapsed os = ITensorMPS.OpSum(Qc)
    mps = mt.get_mps(ψ)
    build = @elapsed W = MPO(os, siteinds(mps))
    inn = @elapsed ev_split = inner(mps', W, mps)
    maxdim_mpo = maxlinkdim(W)
    abs(ev_split - ev_seam) <= 1e-10 || @warn "split run differs from the seam" ev_split ev_seam
    flushprintln("seam $(fmt(seam)) s (OpSum $(fmt(opsum)), MPO $(fmt(build)), inner $(fmt(inn)); MPO maxdim $maxdim_mpo)")
  end

  bft = @benchmark frame_terms($ψ, $Q) samples = 5 evals = 1 seconds = 10
  terms = frame_terms(ψ, Q)
  naive, ev_naive, naive_est = timed_perterm(perterm_naive, terms)
  flushprintln("per-term naive $(fmt(naive)) s (estimate $(fmt(naive_est)) s)")
  prefix, ev_prefix, prefix_est = timed_perterm(perterm_prefix, terms)
  flushprintln("per-term prefix $(fmt(prefix)) s (estimate $(fmt(prefix_est)) s), frame_terms $(fmt(med(bft))) s")

  # Agreement: per-term variants with each other and with the seam; the seam with the old value.
  ref = isnan(real(ev_seam)) ? ev_prefix : ev_seam
  for (name, v) in (("naive", ev_naive), ("prefix", ev_prefix), ("old", ev_old))
    isnan(real(v)) || isnan(real(ref)) || flushprintln(@sprintf("  |ev_%s − ref| = %.1e", name, abs(v - ref)))
  end
  push!(rows, (n, med(bconv), med(bconj), seam, opsum, build, inn, maxdim_mpo, med(bft),
               naive, naive_est, prefix, prefix_est, ev_seam, ev_naive, ev_prefix, ev_old))
end

mkpath(dirname(outpath))
open(outpath, "w") do io
  println(io, "# Phase 5 new-stack benchmarks (bench_new.jl), $(gethostname()), Julia $(VERSION), ",
          "threads $(Threads.nthreads()), BLAS $(BLAS.get_num_threads()); state bond dimension $bd")
  println(io, "# times in seconds: BenchmarkTools medians, except seam/opsum/mpo_build/inner and per-term runs estimated above 60 s (one run each); ",
          "NaN = not run (above maxmpo, or estimate above budget = $budget s)")
  println(io, join(cols, '\t'))
  for r in rows
    println(io, join([x isa Complex ? repr(x) : string(x) for x in r], '\t'))
  end
end
flushprintln("\nwrote $outpath")
