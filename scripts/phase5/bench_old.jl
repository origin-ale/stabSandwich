# Phase 5: local timings of the OLD stack (vendored CliffordMPS), from the `pre-port` worktree's
# environment: the old pp → cmps.PauliSum converter and the old `expectation` (conjugate → MPO).
# Single @elapsed runs, not BenchmarkTools (an old MPO evaluation at 10³ terms takes minutes).
# Run from this repository's root with the worktree's project, under a timeout:
#
#   julia -t 1 --project=../stabSandwich-old scripts/phase5/bench_old.jl [maxold=1161] [budget=300] [mpobudget=600] [out=output/phase5/bench_old.tsv]
#
#   maxold     largest sum timed on the old stack (plan: old timings up to about 10³ terms)
#   budget     seconds allowed for the old term-by-term run per sum (estimated from its first terms)
#   mpobudget  seconds allowed for one old MPO evaluation, estimated from the previous size as
#              t_prev·(n/n_prev)^2.75 (the random-sum scaling; real sums compressed better
#              in T8); a sum over budget is skipped. The first sum (201 terms, 3.8 s in T8) always runs.
#
# Converter: each size is converted twice (first and second call). The plan's hypothesis is that
# the old converter's `vcat(xs...)` compiles a method per term count, so the first call at a new
# size is slow and the second fast. The probe at the end converts a different sum of an already
# seen size and one of a new size, to separate "new size" from "new content".

using CampsPP, CampsPP.DisentangleCAMPS
import CliffordMPS as cmps
using QuantumClifford
using LinearAlgebra, Strided
Strided.disable_threads()
BLAS.set_num_threads(1)
include("common.jl")

args = Dict(split(a, "="; limit = 2) |> kv -> kv[1] => kv[2] for a in ARGS)
for k in keys(args)
  k in ("maxold", "budget", "mpobudget", "out") || error("Unrecognised argument `$k` (expected maxold, budget, mpobudget, out)")
end
maxold = parse(Int, get(args, "maxold", "1161"))
budget = parse(Float64, get(args, "budget", "300"))
mpobudget = parse(Float64, get(args, "mpobudget", "600"))
outpath = get(args, "out", joinpath(@__DIR__, "..", "..", "output", "phase5", "bench_old.tsv"))
flushprintln("project: $(Base.active_project())")
occursin("stabSandwich-old", Base.active_project()) || error("run with --project=<the pre-port worktree>")
flushprintln("threads: julia $(Threads.nthreads()), BLAS $(BLAS.get_num_threads()); maxold = $maxold, budget = $budget s, mpobudget = $mpobudget s")

# == State (as generate.jl: CAMPS untruncated until the bond dimension reaches χ) ==
st = read_state()
N = st.N
gates, phases = sample_circuit(st, xxz_circuit)
O = transferredmagnetization(N, st.excited)
t_state = @elapsed (ψ, k, s, _) = camps_circuit_dynamics(computationalcamps(N, st.excited), gates, phases, CHI, O)
s == st.s || error("rebuild switched at gate $s, expected $(st.s)")
bd = DisentangleCAMPS.bonddim(ψ)
flushprintln("state: s = $s, k = $k, bond dimension $bd ($(round(t_state; digits = 1)) s)")

strop(str) = PauliOperator(0x0, [c in "XY" for c in str], [c in "ZY" for c in str])
termwise(ψ, rows) = sum(c * cmps.expectation(ψ, strop(str)) for (str, c) in rows)

nterms_of(f) = parse(Int, match(r"sum_0*(\d+)\.txt", basename(f))[1])
files = filter(f -> nterms_of(f) <= maxold, sum_files())
# Warm-up on a 3-term sum: compiles the generic code paths, and vcat only for 3 terms.
let rows = read_sum(first(files))[3][1:3]
  S3 = ppsum_from_rows(N, rows)
  cmps.expectation(ψ, cmps.PauliSum(S3)); termwise(ψ, rows)
end

fmt(x) = isfinite(x) ? @sprintf("%.4g", x) : "—"
cols = ["nterms", "conv1_s", "conv2_s", "mpo_s", "mpo_est_s", "termwise_s", "termwise_est_s", "ev_mpo", "ev_termwise", "ev_t8_old"]
results = []
prev_mpo = nothing                             # (n, seconds) of the last old MPO run
for f in files
  S, header, rows = read_sum(f)
  n = length(S)
  flushprintln("\n== $(basename(f)): $n terms")
  GC.gc()
  conv1 = @elapsed op = cmps.PauliSum(S)
  conv2 = @elapsed cmps.PauliSum(S)
  flushprintln("converter: first call $(fmt(conv1)) s, second $(fmt(conv2)) s")
  GC.gc()
  mpo_est = prev_mpo === nothing ? 0.0 : prev_mpo[2] * (n / prev_mpo[1])^2.75
  if mpo_est <= mpobudget
    mpo = @elapsed ev_mpo = cmps.expectation(ψ, op)
    global prev_mpo = (n, mpo)
  else
    mpo, ev_mpo = NaN, NaN + 0im
  end
  flushprintln("expectation (conjugate + MPO): $(fmt(mpo)) s (estimate $(fmt(mpo_est)) s), ⟨O⟩ = $ev_mpo")
  n0 = min(20, n)
  t0 = @elapsed termwise(ψ, rows[1:n0])
  est = t0 * n / n0
  tw, ev_tw = est <= budget ? (@elapsed(v = termwise(ψ, rows)), v) : (NaN, NaN + 0im)
  flushprintln("term by term: $(fmt(tw)) s (estimate $(fmt(est)) s)")
  ev_t8 = haskey(header, "old_expectation_mpo") ?
          parse(ComplexF64, replace(header["old_expectation_mpo"], " " => "")) : NaN + 0im
  isnan(real(ev_t8)) || flushprintln(@sprintf("  |ev_mpo − T8 old| = %.1e", abs(ev_mpo - ev_t8)))
  isnan(real(ev_tw)) || isnan(real(ev_mpo)) || flushprintln(@sprintf("  |ev_termwise − ev_mpo| = %.1e", abs(ev_tw - ev_mpo)))
  push!(results, (n, conv1, conv2, mpo, mpo_est, tw, est, ev_mpo, ev_tw, ev_t8))
end

# == Converter compilation probe ================================================
# Rows of the 538-term sum, cut to sizes already converted (201) and not yet seen (202, 203).
flushprintln("\n== converter probe (rows of sum_0538 cut to a given size)")
rows538 = read_sum(joinpath(T8DIR, "sum_0538.txt"))[3]
probe = []
for n in (201, 202, 202, 203, 538)
  Sn = ppsum_from_rows(N, rows538[1:n])
  t = @elapsed cmps.PauliSum(Sn)
  flushprintln("  $n terms: $(fmt(t)) s")
  push!(probe, (n, t))
end

mkpath(dirname(outpath))
open(outpath, "w") do io
  println(io, "# Phase 5 old-stack timings (bench_old.jl, project $(Base.active_project())), $(gethostname()), ",
          "Julia $(VERSION), threads $(Threads.nthreads()), BLAS $(BLAS.get_num_threads()); state bond dimension $bd")
  println(io, "# single @elapsed runs in seconds; NaN = not run (estimate above budget: term by term $budget s, MPO $mpobudget s)")
  println(io, join(cols, '\t'))
  for r in results
    println(io, join([x isa Complex ? repr(x) : string(x) for x in r], '\t'))
  end
  println(io, "# converter probe (rows of sum_0538 cut to n, in this order): n, seconds")
  for (n, t) in probe
    println(io, "# probe\t$n\t$t")
  end
end
flushprintln("\nwrote $outpath")
