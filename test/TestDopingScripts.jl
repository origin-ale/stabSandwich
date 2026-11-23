# Post-port item 9, script level: count-mode doping in Comptime_prob_tm.jl and
# TransferredMag_varparams.jl (test/CONVENTIONS.md §2 "Doping": count mode, seeds in the
# scripts, output labels). Each script runs small, in child processes on the root
# environment, and its output files are checked against the spec:
#   - Comptime_prob_tm.jl, all four methods, n = 0, 3 and T = 28 (`:xy`), and CAMPS-PP alone
#     in `:z` mode (n = 3, 28): every sample's doped gates are exactly n XX-YY pairs (`:xy`)
#     or n ZZ gates (`:z`) on n distinct triples, and those triples are the mask of
#     MersenneTwister(hash((rng_seed + i, n))); initial states from seed rng_seed + i; the x
#     column of comptime_prob_tm.txt is n; every block of every method is headed
#     `# n = <integer>, μ = <μ>`. What ran is checked too, not only what was recorded: each
#     method's per-sample ⟨O⟩ equals the dense reference on the spec's state and mask.
#   - TransferredMag_varparams.jl: a (μ, sample) has the same initial state for every n
#     (state stream on (μ index, sample) only), states independent across μ (not the same
#     uniforms thresholded at each μ), masks independent across n (not nested, as a shared
#     state/doping stream would make them), across μ and across samples, exactly n triples
#     each, n = T included.
#   - Both scripts refuse n = T + 1 and count lists that are not strictly increasing, before
#     writing any output.
#   - Comptime_prob_tm.jl's per-method caps (campspp_nmax, camps_nmax, mps_nmax, pp_nmax):
#     a method runs its samples at the counts n ≤ its cap only, and at the others has NaN
#     in its seven columns of comptime_prob_tm.txt and gctime_prob_tm.txt, no block in its
#     per-n files and no sample run; a cap below every count is refused before any output.
#   - Comptime_prob_tm.jl's timing summaries leave out the samples that stopped early on
#     Pmax (CAMPS-PP, PP): each method's seven columns of comptime_prob_tm.txt and
#     gctime_prob_tm.txt are the mean, std. err. (over the kept count), quartiles and extrema
#     of the per-sample times of its full-trajectory samples, all NaN if none is; CAMPS and
#     MPS keep every sample; the per-sample timing files still list every sample; a warning
#     names exactly each stopped sample.
# The spec fixes the Comptime_prob_tm.jl seeds exactly; for TransferredMag_varparams.jl it
# fixes only what each stream depends on, so those checks are relational. Each relational
# check asserts that the signature of the wrong seeding (e.g. "nested for every sample")
# does not hold for every sample; the comment at each check bounds the probability that
# independent draws show the signature everywhere, so none of them passes by luck.
#
#   julia --project=. -t 2 test/TestDopingScripts.jl   (also run by test/runtests.jl)
# About 7 min: five script runs of 30–110 s (mostly compilation) and seven refused runs of 5 s.
#
# STABSANDWICH_SCRIPT_DIR (default: the repository root) points the runs at another copy of
# the scripts; used for the scratch mutants of item 9, not by the suite.

using Test
using CampsPP
using CampsPP.DisentangleCAMPS
using Random: MersenneTwister
using Statistics: mean, std, quantile

const ROOT = dirname(@__DIR__)
const SCRIPT_DIR = get(ENV, "STABSANDWICH_SCRIPT_DIR", ROOT)

isdefined(@__MODULE__, :DenseReference) ||
  include(joinpath(ROOT, "CampsPP", "test", "DisentangleCAMPS", "DenseReference.jl"))
const DR = DenseReference

"Run `script` in SCRIPT_DIR with `args` in a child process; returns (exit code, stdout, stderr)."
function run_script(script, args...; threads = 1)
  cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no -t $threads $(joinpath(SCRIPT_DIR, script)) $args`
  obuf, ebuf = IOBuffer(), IOBuffer()
  proc = run(pipeline(ignorestatus(cmd); stdout = obuf, stderr = ebuf))
  return proc.exitcode, String(take!(obuf)), String(take!(ebuf))
end

const LABEL_RE = r"^# n = (\d+), μ = (\S+?)(, max sample = \d+)?$"

"""
Split `path` into its labelled blocks: a vector of `(n, μ, lines)`, one per `# n = …, μ = …`
header, with the block's other non-empty lines. Lines before the first header are dropped.
"""
function label_blocks(path)
  blocks = Tuple{Int, String, Vector{String}}[]
  for l in eachline(path)
    m = match(LABEL_RE, l)
    if m !== nothing
      push!(blocks, (parse(Int, m[1]), m[2], String[]))
    elseif !isempty(blocks) && !isempty(strip(l))
      push!(blocks[end][3], l)
    end
  end
  return blocks
end

"Lines of `path` that look like a block label (`# <word> = …`)."
label_lines(path) = filter(l -> occursin(r"^# \w+ = ", l), readlines(path))

"""
Doped gates of `doped_prob_tm.txt`: `Dict(n => [gates of sample 1, gates of sample 2, …])`,
each gate `(index, Pauli string)`; also checks each `# sample s: k doped gates` count line
against the gates listed under it (returns `nothing` on a mismatch).
"""
function doped_gates(path)
  out = Dict{Int, Vector{Vector{Tuple{Int, String}}}}()
  for (n, _, lines) in label_blocks(path)
    samples = Vector{Tuple{Int, String}}[]
    declared = Int[]
    for l in lines
      m = match(r"^# sample (\d+): (\d+) doped gates$", l)
      if m !== nothing
        parse(Int, m[1]) == length(samples) + 1 || return nothing
        push!(samples, Tuple{Int, String}[]); push!(declared, parse(Int, m[2]))
      else
        f = split(l)
        push!(samples[end], (parse(Int, f[1]), String(f[3])))
      end
    end
    length.(samples) == declared || return nothing
    out[n] = samples
  end
  return out
end

"""
Triples doped by the gate list `gs`, or `nothing` unless the gates are exactly an XX at 3m−2
and a YY at 3m−1 for each listed triple m (`:xy`), or exactly a ZZ at 3m (`:z`).
"""
function doped_triples(gs, mode)
  ms = sort(unique([cld(g, 3) for (g, _) in gs]))
  want = mode == :xy ? sort([[(3m - 2, "XX") for m in ms]; [(3m - 1, "YY") for m in ms]]) :
                       [(3m, "ZZ") for m in ms]
  sort([(g, filter(!=('I'), P)) for (g, P) in gs]) == want || return nothing
  return ms
end

"""
Per-sample ⟨O⟩ of an `evs_prob_tm_<method>_<thl>_full.txt` file:
`Dict(n => [values of sample 1, values of sample 2, …])`, layer 0 first.
"""
function sample_evs(path)
  out = Dict{Int, Vector{Vector{Float64}}}()
  for (n, _, lines) in label_blocks(path)
    samples = Vector{Float64}[]
    for l in lines
      if startswith(l, "# sample ")
        push!(samples, Float64[])
      else
        push!(samples[end], parse(Float64, split(l)[2]))
      end
    end
    out[n] = samples
  end
  return out
end

# -- Dense reference for the scripts' circuit (CONVENTIONS.md §1, §2, §5) -----------------
# Built from the spec alone: bricklayer layers (bonds (1,2), (3,4), … then (2,3), …), each
# bond an XX, YY, ZZ triple in the exp(iϕP) convention; ⟨O⟩ of the transferred
# magnetization O = −½ Σ_L Z_j + ½ Σ_R Z_j + (N/2 − n₁ᴸ(0) − n₀ᴿ(0)), recorded at t = 0 and at
# each layer end 3(N−1)·l.

"Gate strings of t bricklayer XXZ layers on N qubits, in gate order (§2)."
function spec_gate_strings(N, t)
  bonds = vcat([(i, i + 1) for i in 1:2:N-1], [(i, i + 1) for i in 2:2:N-1])
  pair(P, (a, b)) = String([q in (a, b) ? P : 'I' for q in 1:N])
  return repeat([pair(P, b) for b in bonds for P in ('X', 'Y', 'Z')], t)
end

"""
Angles of the doped circuit: ϕ_xy on XX, YY and ϕ_z on ZZ, then the doped triples `ms` get
`magicphase` on XX and YY (`:xy`) or on ZZ (`:z`), as `xy_magic`/`z_magic` do (§2).
"""
function spec_phases(N, t, ms, mode; ϕ_xy, ϕ_z, magicphase)
  ph = repeat([ϕ_xy, ϕ_xy, ϕ_z], (N - 1) * t)
  for m in ms
    mode == :xy ? (ph[3m - 2] = ph[3m - 1] = magicphase) : (ph[3m] = magicphase)
  end
  return ph
end

"Dense ⟨O⟩ at t = 0 and at each layer end, from the basis state with `ones` in |1⟩ (§5)."
function dense_tm(N, t, ones, phases)
  L, R = 1:N÷2, N÷2+1:N
  z(q) = String([j == q ? 'Z' : 'I' for j in 1:N])
  c0 = N / 2 - count(in(ones), L) - count(!in(ones), R)
  O(ψ) = -sum(DR.expect(ψ, z(q)) for q in L) / 2 + sum(DR.expect(ψ, z(q)) for q in R) / 2 + c0
  ψ = DR.basis_state(N, ones)
  evs = [O(ψ)]
  for (k, (P, ϕ)) in enumerate(zip(spec_gate_strings(N, t), phases))
    DR.rotate!(ψ, P, ϕ)
    k % (3 * (N - 1)) == 0 && push!(evs, O(ψ))
  end
  return evs
end

maxabsdiff(a, b) = length(a) == length(b) ? maximum(abs.(a .- b)) : Inf

@testset "dense helpers (hand pins)" begin
  # N=4, t=1: bonds (1,2), (3,4), (2,3); one ZZ on (2,3) is the 9th gate
  @test spec_gate_strings(4, 1)[[1, 2, 3, 4, 9]] == ["XXII", "YYII", "ZZII", "IIXX", "IZZI"]
  @test spec_phases(4, 2, [2], :xy; ϕ_xy = 0.1, ϕ_z = 0.2, magicphase = 0.5)[[4, 5, 6]] == [0.5, 0.5, 0.2]
  @test spec_phases(4, 2, [2], :z; ϕ_xy = 0.1, ϕ_z = 0.2, magicphase = 0.5)[[4, 5, 6]] == [0.1, 0.1, 0.5]
  # ⟨O⟩ = 0 at t = 0 (§5), for any basis state
  @test all(dense_tm(4, 0, ones, Float64[])[1] == 0 for ones in ([], [1], [3, 4], [1, 2, 3]))
  # XXZ at ϕ_xy = ϕ_z = π/4 is a SWAP (§2, TestCircuits.jl): on |1000⟩ (N=4, t=1) bonds
  # (1,2), (3,4), (2,3) move the excitation 1 → 2 → 3, from L to R: ⟨O⟩ = −2 (§5)
  @test dense_tm(4, 1, [1], fill(π/4, 9)) ≈ [0, -2] atol = 1e-12
end

# -- Comptime_prob_tm.jl ------------------------------------------------------------------
# All four methods at N = 8 (T = (N−1)·N/2 = 28 triples, N/2 = 4 layers), in the exact regime,
# so each method's ⟨O⟩ is compared with the dense reference at 1e-10 (§6):
#   - CAMPS-PP: χ = 64 is never reached (bond dimension ≤ 2^(N/2) = 16), so it never switches
#     to PP and its thl is unused; checked through campspp_switch_prob_tm.txt (84 = never);
#   - MPS: mps_thl = 0, no truncation (it renormalises a state of norm 1);
#   - PP: pp_thl = 0 and pp_Pmax ≥ 4^8, no truncation;
#   - CAMPS: camps_thl = 1e-11, the smallest threshold §2 allows for tests that depend on the
#     layer-end disentangler (below 1e-13 the flag is unreliable). This is nominally the
#     truncated regime (§6 gives no bound there); see CT_ATOL.
const CT_N, CT_SAMPLES, CT_SEED, CT_μ = 8, 2, 100, 0.6
const CT_T = (CT_N - 1) * (CT_N ÷ 2)
const CT_LAYERS = CT_N ÷ 2
const CT_METHODS = ["campspp", "camps", "mps", "pp"]
const CT_ARGS = ["N=$CT_N", "samples=$CT_SAMPLES", "methods=$(join(CT_METHODS, ','))",
                 "camps_thl=1e-11", "mps_thl=0", "pp_thl=0"]
# Script defaults (checked against the header): ϕ = θ = π/4; magic angle 3π/16 on XX-YY
# (`:xy`), π/3 on ZZ (`:z`), the values the script's parameter comment states.
const CT_MAGIC = Dict(:xy => 3π/16, :z => π/3)
# Exact regime, §6: 1e-10 for every method. For CAMPS at camps_thl = 1e-11 (nominally
# truncated) the measured max |⟨O⟩ − dense| over this run is 1.4e-14 (2026-11-20), so the
# exact-regime tolerance is kept for it too, not loosened.
const CT_ATOL = 1e-10

# §2: sample s's mask at count n is doping_mask_count(MersenneTwister(hash((rng_seed + s, n))), …),
# its initial state from seed rng_seed + s, the same for every n (s = 0: the warm-up sample)
ct_spec_mask(n, s) = findall(doping_mask_count(MersenneTwister(hash((CT_SEED + s, n))),
                                               last(xxz_circuit(π/4, π/4, CT_LAYERS, CT_N)), n))
ct_spec_ones(s) = domainwallstate(MersenneTwister(CT_SEED + s), CT_N, CT_μ)[2]
"Dense ⟨O⟩ of sample `s` at count `n` in magic mode `mode`, at the script defaults."
ct_dense(n, s, mode) = dense_tm(CT_N, CT_LAYERS, ct_spec_ones(s),
  spec_phases(CT_N, CT_LAYERS, ct_spec_mask(n, s), mode; ϕ_xy = π/4, ϕ_z = π/4,
              magicphase = CT_MAGIC[mode]))

"""
Spec checks on a Comptime_prob_tm.jl output directory run with `counts`, `methods` and
magic mode `mode`; returns a Dict of named Bools.
"""
function check_comptime(dir; counts, methods, mode)
  N, t = CT_N, CT_LAYERS
  T = CT_T
  r = Dict{String, Bool}()
  header = readline(joinpath(dir, "comptime_prob_tm.txt"))
  r["header lists magic_counts, ϕ, θ, magic_mode"] =
    occursin("magic_counts=$(join(counts, ','))", header) &&
    occursin("ϕ=$(π/4) θ=$(π/4)", header) && occursin("magic_mode=$mode", header)

  spec_mask, spec_ones = ct_spec_mask, ct_spec_ones

  doped = doped_gates(joinpath(dir, "doped_prob_tm.txt"))
  r["doped file parses"] = doped !== nothing && sort(collect(keys(doped))) == counts &&
                           all(length(doped[n]) == CT_SAMPLES for n in counts)
  if r["doped file parses"]
    ms = Dict((n, s) => doped_triples(doped[n][s], mode) for n in counts, s in 1:CT_SAMPLES)
    r["doped gates: n triples, $(mode == :xy ? "2n XX-YY" : "n ZZ") gates"] =
      all(v !== nothing && length(v) == n && all(1 .<= v .<= T) for ((n, _), v) in ms)
    r["mask from hash((rng_seed + i, n))"] =
      all(ms[(n, s)] == spec_mask(n, s) for n in counts, s in 1:CT_SAMPLES)
  end

  # §2: sample i's initial state from seed rng_seed + i, the same for every n
  states = [split(l) for l in readlines(joinpath(dir, "initstates_prob_tm.txt"))
            if !isempty(l) && !startswith(l, '#')]
  r["initial states from rng_seed + i"] = length(states) == CT_SAMPLES && all(
    parse(Int, f[1]) == s && parse.(Int, f[3:end]) == spec_ones(s) for (s, f) in enumerate(states))

  # §2: the x column of comptime_prob_tm.txt (and gctime) is n, written as an integer
  r["x column is n"] = all(begin
      rows = [split(l) for l in readlines(joinpath(dir, f)) if !isempty(l) && !startswith(l, '#')]
      [r[1] for r in rows] == string.(counts)
    end for f in ("comptime_prob_tm.txt", "gctime_prob_tm.txt"))

  # §2: every block label in every output, of every method, is `# n = <integer>, μ = <μ>`,
  # for each n in order
  files = filter(f -> endswith(f, ".txt") && !endswith(f, "_log.txt"), readdir(dir))
  labels = [(f, l) for f in files for l in label_lines(joinpath(dir, f))]
  r["labels present for every method"] = all(
    any(occursin("_$m", f) for (f, _) in labels) for m in methods)
  r["every label is # n = <int>, μ = 0.6"] = all(occursin(LABEL_RE, l) &&
    match(LABEL_RE, l)[2] == string(CT_μ) for (_, l) in labels)
  r["blocks in n order"] = all(begin
      ns = [b[1] for b in label_blocks(joinpath(dir, f))]
      isempty(ns) || unique(ns) == counts
    end for f in files)

  # What ran, not only what was recorded (review follow-ups 4 and 5): each method's ⟨O⟩ of
  # sample s at count n equals the dense ⟨O⟩ on the spec's state and mask for (s, n). This is
  # also "the same mask for every method", block by block.
  dense = Dict((n, s) => ct_dense(n, s, mode)
               for n in counts, s in 0:CT_SAMPLES)          # s = 0: the warm-up sample
  for m in methods
    evs = sample_evs(joinpath(dir, "evs_prob_tm_$(m)_10_full.txt"))
    errs = [maxabsdiff(get(evs, n, [])[s], dense[(n, s)])
            for n in counts, s in 1:CT_SAMPLES if length(get(evs, n, [])) == CT_SAMPLES]
    ok = length(errs) == length(counts) * CT_SAMPLES && maximum(errs) <= CT_ATOL
    ok || println("  $m: max |⟨O⟩ − dense| = ", isempty(errs) ? "no data" : maximum(errs))
    m == "camps" && println("  camps (camps_thl = 1e-11): max |⟨O⟩ − dense| = ",
                            isempty(errs) ? "no data" : maximum(errs))
    r["$m: ⟨O⟩ = dense on the spec's state and mask"] = ok
  end
  # Negative check: the comparison tells samples apart. For every n, some sample's dense
  # ⟨O⟩ differs by > 0.1 from that of the sample before it (the mutant that runs sample
  # i − 1's state and mask), and in `:xy` mode, for every doped n, some sample's differs by
  # > 0.1 from the same state with no doping (a method that ignores the mask). Not in `:z`
  # mode: at ϕ_xy = π/4 every XX-YY pair maps basis states to basis states, so ZZ gates
  # only add phases and ⟨O⟩ does not depend on a `:z` mask. There the mask is checked
  # through the doped gates only.
  r["negative: dense(s) differs from dense(s − 1)"] = all(
    any(maxabsdiff(dense[(n, s)], dense[(n, s - 1)]) > 0.1 for s in 1:CT_SAMPLES) for n in counts)
  if mode == :xy
    r["negative: dense(s) differs without the mask"] = all(
      any(maxabsdiff(dense[(n, s)], dense_tm(N, t, spec_ones(s), fill(π/4, 3T))) > 0.1
          for s in 1:CT_SAMPLES) for n in counts if n > 0)
  end

  # CAMPS-PP never switched to PP (column 2 = 84 = number of gates: never), so it ran in
  # the exact regime
  if "campspp" in methods
    rows = [split(l) for l in readlines(joinpath(dir, "campspp_switch_prob_tm.txt"))
            if !isempty(strip(l)) && !startswith(l, '#')]
    r["campspp never switched (exact regime)"] =
      length(rows) == length(counts) * CT_SAMPLES && all(f[2] == string(3T) for f in rows)
  end
  return r
end

# -- Comptime_prob_tm.jl, per-method caps --------------------------------------------------
# One run at the CT_ARGS settings (exact regime) with counts 0, 3, 10, T, and a different
# run set per method, so a cap read from another method's parameter, or compared with < for
# ≤, shows up:
#   campspp: its default, the largest count (of the magic_counts override, not of the
#            script's own list) → every n;
#   camps:   10, a count (n = cap runs) → 0, 3, 10;
#   mps:     0, the first count → 0 only;
#   pp:      5, between counts → 0, 3.
const NMAX_COUNTS = [0, 3, 10, CT_T]
const NMAX_CAPS = Dict("camps" => 10, "mps" => 0, "pp" => 5)    # campspp: the default
const NMAX_PROGRESS = Dict("CAMPS-PP" => "campspp", "CAMPS" => "camps", "MPS" => "mps",
                           "PP" => "pp")

"""
Checks on a Comptime_prob_tm.jl output directory and stdout from the run with NMAX_COUNTS
and NMAX_CAPS; returns a Dict of named Bools.
"""
function check_nmax(dir, stdout_)
  counts = NMAX_COUNTS
  caps = merge(Dict("campspp" => maximum(counts)), NMAX_CAPS)
  runs = Dict(m => filter(<=(caps[m]), counts) for m in CT_METHODS)
  r = Dict{String, Bool}()

  header = readline(joinpath(dir, "comptime_prob_tm.txt"))
  r["header lists each cap (campspp: the largest count)"] =
    all(occursin(Regex("\\b$(m)_nmax=$(caps[m])\\b"), header) for m in CT_METHODS)

  # Seven columns per method, in `methods` order after the n column: all finite where the
  # method ran, all NaN where it did not
  for f in ("comptime_prob_tm.txt", "gctime_prob_tm.txt")
    rows = [parse.(Float64, split(l)) for l in readlines(joinpath(dir, f))
            if !isempty(l) && !startswith(l, '#')]
    shape_ok = [row[1] for row in rows] == counts &&
               all(length(row) == 1 + 7 * length(CT_METHODS) for row in rows)
    r["$f: one row per n, 1 + 7·4 fields"] = shape_ok
    shape_ok || continue
    for (k, m) in enumerate(CT_METHODS)
      cols = 1 + 7 * (k - 1) .+ (1:7)
      r["$f, $m: finite at n ≤ $(caps[m]), NaN above"] = all(
        n in runs[m] ? all(isfinite, row[cols]) : all(isnan, row[cols])
        for (n, row) in zip(counts, rows))
    end
  end

  # Every labelled file of method m (name has `_<m>_` or `_<m>.`, plus CAMPS-PP's switch
  # file) has blocks at exactly the counts it ran; the doped gates are still written for
  # every n
  files = filter(f -> endswith(f, ".txt") && !endswith(f, "_log.txt"), readdir(dir))
  for m in CT_METHODS
    mfiles = filter(f -> occursin(Regex("_$(m)[_.]"), f) ||
                         (m == "campspp" && f == "campspp_switch_prob_tm.txt"), files)
    r["$m: per-n files found"] = length(mfiles) >= 3    # MPS: times, evs, evs full
    r["$m: blocks exactly at n ≤ $(caps[m])"] = all(
      unique(b[1] for b in label_blocks(joinpath(dir, f))) == runs[m] for f in mfiles)
  end
  doped = doped_gates(joinpath(dir, "doped_prob_tm.txt"))
  r["doped gates written for every n"] = doped !== nothing && sort(collect(keys(doped))) == counts

  # What ran is still right (skipping a method does not shift the others' masks or states):
  # per-sample ⟨O⟩ = dense at each count the method ran
  for m in CT_METHODS
    evs = sample_evs(joinpath(dir, "evs_prob_tm_$(m)_10_full.txt"))
    r["$m: ⟨O⟩ = dense at each n it ran"] = sort(collect(keys(evs))) == runs[m] && all(
      length(evs[n]) == CT_SAMPLES &&
      all(maxabsdiff(evs[n][s], ct_dense(n, s, :xy)) <= CT_ATOL for s in 1:CT_SAMPLES)
      for n in runs[m])
  end

  # The work itself is skipped, not only its output: the progress log (frac = 0, one line
  # per sample, on stdout) lists samples 1…samples of m at exactly the counts it ran, and
  # each skip is logged
  ran = Dict(m => Dict{Int, Vector{Int}}() for m in CT_METHODS)
  for mt in eachmatch(r"\] n=(\d+) (CAMPS-PP|CAMPS|MPS|PP) sample (\d+)/", stdout_)
    push!(get!(ran[NMAX_PROGRESS[mt[2]]], parse(Int, mt[1]), Int[]), parse(Int, mt[3]))
  end
  for m in CT_METHODS
    r["$m: samples run at exactly n ≤ $(caps[m])"] = sort(collect(keys(ran[m]))) == runs[m] &&
      all(ran[m][n] == 1:CT_SAMPLES for n in runs[m])
    r["$m: skip logged at each n > $(caps[m])"] = all(
      occursin("n=$n: skipping $m (above $(m)_nmax = $(caps[m]))", stdout_)
      for n in setdiff(counts, runs[m]))
  end
  return r
end

# -- Comptime_prob_tm.jl, samples stopped on Pmax -------------------------------------------
# One run at N = 8 with 6 samples, counts 0, 3, 6, T, all four methods, every truncation
# threshold 0 except CAMPS' (as in CT_ARGS), so which samples stop is deterministic:
#   campspp: χ = 2 makes it switch to PP at a sample-dependent gate, and Pmax = 2000 then
#            stops some samples: none at n = 0 (no doping, never switches), some at n = 3
#            and 6, all at n = T (2026-11-23: 0, 1, 1, 6 stopped);
#   pp:      Pmax = 2000 stops every sample at every n (at pp_thl = 0 the Pauli count does
#            not depend on the state or mask: 9929 terms after layer 2 of 4);
#   camps:   capped at n = 3, so it runs at a count where CAMPS-PP has both kinds of sample
#            (a block that took CAMPS-PP's kept flags would show there) and stays short;
#   mps:     every n.
# Which samples stopped is read from what ran (an ⟨O⟩ trajectory short of N/2 + 1 points),
# not pinned, and the scenario's coverage is checked so the summary checks are not vacuous.
const PMAX_SAMPLES = 6
const PMAX_COUNTS = [0, 3, 6, CT_T]
const PMAX_CAPS = Dict("campspp" => CT_T, "camps" => 3, "mps" => CT_T, "pp" => CT_T)
const PMAX_ARGS = ["N=$CT_N", "samples=$PMAX_SAMPLES", "methods=$(join(CT_METHODS, ','))",
                   "magic_counts=$(join(PMAX_COUNTS, ','))", "camps_nmax=3",
                   "camps_thl=1e-11", "mps_thl=0", "pp_thl=0", "campspp_thl=0",
                   "campspp_χ=2", "campspp_Pmax=2000", "pp_Pmax=2000"]
const PMAX_PROGRESS = Dict("CAMPS-PP" => "campspp", "PP" => "pp")

"""
Per-sample times of a `comptime_full_prob_tm_<method>.txt` file: `Dict(n => rows)`, each row
`[sample, time, GC time]`.
"""
sample_times(path) = Dict(n => [parse.(Float64, split(l)) for l in lines]
                          for (n, _, lines) in label_blocks(path))

"""
The seven summary columns of `xs` (§ Comptime_prob_tm.jl output): mean, std. err. (std over
√length(xs)), q1, median, q3, min, max; all NaN if `xs` is empty.
"""
time_summary(xs) = isempty(xs) ? fill(NaN, 7) :
  [mean(xs), std(xs) / sqrt(length(xs)), quantile(xs, [0.25, 0.5, 0.75])...,
   minimum(xs), maximum(xs)]

"Elementwise equal up to rounding (values are written at full precision), NaN = NaN."
same_summary(a, b) = length(a) == length(b) &&
  all(isnan(x) ? isnan(y) : isapprox(x, y; rtol = 1e-12) for (x, y) in zip(a, b))

"""
Checks on a Comptime_prob_tm.jl output directory and stdout from the run with PMAX_ARGS;
returns a Dict of named Bools.
"""
function check_pmax(dir, stdout_)
  S, counts = PMAX_SAMPLES, PMAX_COUNTS
  runs = Dict(m => filter(<=(PMAX_CAPS[m]), counts) for m in CT_METHODS)
  r = Dict{String, Bool}()

  # What stopped, from what ran: a sample whose ⟨O⟩ trajectory misses the final layer
  stopped = Dict{Tuple{String, Int}, Vector{Int}}()
  times = Dict{String, Dict{Int, Vector{Vector{Float64}}}}()
  for m in CT_METHODS
    evs = sample_evs(joinpath(dir, "evs_prob_tm_$(m)_exact_full.txt"))
    times[m] = sample_times(joinpath(dir, "comptime_full_prob_tm_$m.txt"))
    r["$m: ⟨O⟩ of every sample at each n it ran"] =
      sort(collect(keys(evs))) == runs[m] && all(length(evs[n]) == S for n in runs[m])
    # Per-sample timing files are unchanged: every sample, the stopped ones included
    r["$m: per-sample times list every sample at each n it ran"] =
      sort(collect(keys(times[m]))) == runs[m] &&
      all([row[1] for row in times[m][n]] == 1:S && all(length.(times[m][n]) .== 3)
          for n in runs[m])
    r["$m: ⟨O⟩ of every sample at each n it ran"] &&
      r["$m: per-sample times list every sample at each n it ran"] || return r
    for n in runs[m]
      stopped[(m, n)] = [s for s in 1:S if length(evs[n][s]) < CT_LAYERS + 1]
    end
  end
  nstopped(m) = [length(stopped[(m, n)]) for n in runs[m]]

  # Scenario coverage: CAMPS-PP has an n with no sample stopped, one with some stopped and
  # at least 2 kept (so the std. err. is finite and its √(kept count) differs from √S), and
  # one with all stopped; PP stops everywhere; CAMPS and MPS never stop
  r["scenario: campspp has n with none, some (≥ 2 kept) and all stopped"] =
    0 in nstopped("campspp") && any(1 <= k <= S - 2 for k in nstopped("campspp")) &&
    S in nstopped("campspp")
  r["scenario: pp stops every sample at every n"] = all(==(S), nstopped("pp"))
  r["scenario: camps and mps never stop"] = all(==(0), [nstopped("camps"); nstopped("mps")])

  # The summaries: per method, the seven columns at each n are the summary of the per-sample
  # times (column 2) or GC times (column 3) of the samples that did not stop; NaN where the
  # method did not run or every sample stopped
  for (f, col) in (("comptime_prob_tm.txt", 2), ("gctime_prob_tm.txt", 3))
    rows = [parse.(Float64, split(l)) for l in readlines(joinpath(dir, f))
            if !isempty(l) && !startswith(l, '#')]
    shape_ok = [row[1] for row in rows] == counts &&
               all(length(row) == 1 + 7 * length(CT_METHODS) for row in rows)
    r["$f: one row per n, 1 + 7·4 fields"] = shape_ok
    shape_ok || continue
    for (k, m) in enumerate(CT_METHODS)
      cols = 1 + 7 * (k - 1) .+ (1:7)
      r["$f, $m: summary over the samples that did not stop"] = all(
        same_summary(row[cols], n in runs[m] ?
          time_summary([times[m][n][s][col] for s in setdiff(1:S, stopped[(m, n)])]) :
          fill(NaN, 7))
        for (n, row) in zip(counts, rows))
    end
  end

  # Negative check: at each CAMPS-PP count with both kinds of sample, the summary over every
  # sample differs from the one over the kept samples in its mean and std. err. (and the
  # kept samples' std over √S differs from it over √(kept count)), so the check above tells
  # the exclusion apart from no exclusion
  mixed = [n for n in counts if 1 <= length(stopped[("campspp", n)]) <= S - 2]
  r["negative: kept-only summary differs from all-sample one (mean, std. err.)"] = all(begin
      ts = [row[2] for row in times["campspp"][n]]
      kept = ts[setdiff(1:S, stopped[("campspp", n)])]
      all_, kept_ = time_summary(ts), time_summary(kept)
      !isapprox(all_[1], kept_[1]; rtol = 1e-6) && !isapprox(all_[2], kept_[2]; rtol = 1e-6) &&
        !isapprox(std(kept) / sqrt(S), kept_[2]; rtol = 1e-6)
    end for n in mixed)

  # A warning names exactly each stopped sample (stdout)
  warned = Dict{Tuple{String, Int}, Vector{Int}}()
  for mt in eachmatch(r"WARNING: (CAMPS-PP|PP) sample (\d+) for n=(\d+) stopped early \(excluded from the timing statistics\)!", stdout_)
    push!(get!(warned, (PMAX_PROGRESS[mt[1]], parse(Int, mt[3])), Int[]), parse(Int, mt[2]))
  end
  r["stdout: a warning names exactly each stopped sample"] =
    Dict(k => sort(v) for (k, v) in warned) == filter(kv -> !isempty(kv[2]), stopped)
  return r
end

# -- TransferredMag_varparams.jl -----------------------------------------------------------
# N = 8, t = N ÷ 2 = 4 layers: T = 28 triples. 20 samples, so that the relational checks
# below have negligible false-failure probabilities (bounds at each check).
const TM_N, TM_SAMPLES = 8, 20
const TM_μS = ["0.3", "0.6"]
const TM_COUNTS = [0, 5, 10, 28]

"Domain-wall probability p(μ) = e^μ / (e^μ + e^−μ) (§1)."
dw_p(μ) = exp(μ) / (exp(μ) + exp(-μ))

"""
Do the excited sets `a` (at μ₁) and `b` (at μ₂ > μ₁) of one sample look like the same
uniforms thresholded at each μ? §1: a left-half qubit is |1⟩ with probability 1 − p and a
right-half one with probability p, p increasing in μ. With shared uniforms, a left qubit
excited at μ₂ is excited at μ₁ too, and a right qubit excited at μ₁ is excited at μ₂ too
(whichever way the threshold is written: u > p or u < 1 − p for the left half), so
L(μ₂) ⊆ L(μ₁) and R(μ₁) ⊆ R(μ₂).
"""
function shared_uniform_nesting(a, b, N)
  L = 1:N÷2
  return issubset(filter(in(L), b), filter(in(L), a)) &&
         issubset(filter(!in(L), a), filter(!in(L), b))
end

@testset "shared_uniform_nesting (negative check)" begin
  # Same seed at both μ (the mutant's shared stream): nested for every seed
  N = TM_N
  @test all(shared_uniform_nesting(domainwallstate(MersenneTwister(k), N, 0.3)[2],
                                   domainwallstate(MersenneTwister(k), N, 0.6)[2], N) for k in 1:200)
  # Independent seeds: not nested for some seed
  @test !all(shared_uniform_nesting(domainwallstate(MersenneTwister(k), N, 0.3)[2],
                                    domainwallstate(MersenneTwister(k + 1000), N, 0.6)[2], N)
             for k in 1:TM_SAMPLES)
  @test shared_uniform_nesting([1, 2, 6], [1, 6, 7], 8)        # L {1,2} ⊇ {1}, R {6} ⊆ {6,7}
  @test !shared_uniform_nesting([1, 6], [1, 2, 6], 8)          # L grew
  @test !shared_uniform_nesting([1, 6, 7], [1, 6], 8)          # R shrank
end

"Spec checks on a TransferredMag_varparams.jl output directory; returns a Dict of named Bools."
function check_transferredmagnetization(dir)
  N, S = TM_N, TM_SAMPLES
  T = (N - 1) * (N ÷ 2)                                    # 28
  stem = joinpath(dir, "TMDxy_$(N)_$(S)")
  r = Dict{String, Bool}()

  # (n, μ) blocks, n fastest, in every block-structured output
  pairs_want = [(n, μ) for μ in TM_μS for n in TM_COUNTS]
  for suffix in ("", "_full", "_states", "_doped")
    r["blocks (n, μ), n fastest: $suffix"] =
      [(b[1], b[2]) for b in label_blocks(stem * suffix * ".txt")] == pairs_want
  end
  r["every label is # n = <int>, μ = …"] = all(
    occursin(LABEL_RE, l) for suffix in ("", "_full", "_states", "_doped")
    for l in label_lines(stem * suffix * ".txt"))

  # Initial states: sample → bitstring and excited qubits, per (n, μ)
  states = Dict((b[1], b[2]) => [split(l) for l in b[3]] for b in label_blocks(stem * "_states.txt"))
  r["states parse"] = all(haskey(states, k) && length(states[k]) == S for k in pairs_want) &&
    all(join(q in parse.(Int, f[3:end]) ? '1' : '0' for q in 1:N) == f[2] &&
        parse(Int, f[1]) == s for v in values(states) for (s, f) in enumerate(v))
  if r["states parse"]
    bits(n, μ) = [f[2] for f in states[(n, μ)]]
    exc(μ, s) = parse.(Int, states[(0, μ)][s][3:end])
    # §2: the state stream depends on (μ index, sample) only, so not on n
    r["state identical across n for each (μ, sample)"] =
      all(bits(n, μ) == bits(TM_COUNTS[1], μ) for n in TM_COUNTS, μ in TM_μS)
    # The state depends on the sample (review follow-up 9: replaces "all S states distinct",
    # which fails for about 8 % of seedings at N = 8, μ = 0.6). Two independent samples are
    # equal with probability (p² + (1 − p)²)^N: 0.030 at μ = 0.6 and 0.0075 at μ = 0.3 (N = 8).
    # So "all S = 20 samples equal" has probability ≤ 0.030^19 < 1e-28 under the spec, and
    # probability 1 for a state stream that ignores the sample.
    r["states not all equal across samples"] = all(!allequal(bits(0, μ)) for μ in TM_μS)
    # §2 (review follow-up 7): states at different μ come from independent draws, not from
    # the same uniforms thresholded at each p(μ). Shared uniforms nest every sample (see
    # shared_uniform_nesting). Under independent draws a qubit breaks the nesting with
    # probability v = p₁(1 − p₂) (left: |1⟩ at μ₂, |0⟩ at μ₁; right: |1⟩ at μ₁, |0⟩ at μ₂),
    # independently over the N qubits: p₁ = p(0.3) = 0.6457, p₂ = p(0.6) = 0.7685, v = 0.1495,
    # so a sample is nested with probability (1 − v)^8 = 0.274, and all 20 with
    # 0.274^20 ≈ 6e-12.
    μ1, μ2 = TM_μS
    @assert parse(Float64, μ1) < parse(Float64, μ2)
    r["states not nested across μ for every sample (independent draws)"] =
      !all(shared_uniform_nesting(exc(μ1, s), exc(μ2, s), N) for s in 1:S)
  end

  # Doped triples: sample → set of triples, per (n, μ)
  doped = Dict((b[1], b[2]) => [parse.(Int, split(l)) for l in b[3]]
               for b in label_blocks(stem * "_doped.txt"))
  r["doped parse"] = all(haskey(doped, k) && length(doped[k]) == S &&
                         [v[1] for v in doped[k]] == 1:S for k in pairs_want)
  if r["doped parse"]
    tr(n, μ, s) = doped[(n, μ)][s][2:end]
    r["exactly n distinct triples in 1:T (n = T included)"] = all(
      length(tr(n, μ, s)) == n && allunique(tr(n, μ, s)) && all(1 .<= tr(n, μ, s) .<= T)
      for n in TM_COUNTS, μ in TM_μS, s in 1:S)
    # §2: the mask stream is keyed on (n, μ index, sample). For independent draws at T = 28,
    # a 5-subset lies inside a 10-subset with probability binomial(10, 5)/binomial(28, 5)
    # = 252/98280 = 2.6e-3, and two n-subsets (n = 5, 10) are equal with probability
    # ≤ 1/binomial(28, 5) ≈ 1e-5. Each check fails only if the wrong-seeding signature holds
    # for every one of 2·20 (μ, sample) pairs, or all 20 samples of a (n, μ): ≤ (2.6e-3)^40,
    # (1e-5)^20 and (1e-5)^19 under the spec.
    #   - shared state/doping stream (masks are prefixes of one draw): nested for every pair
    r["masks not nested across n (for every sample)"] = !all(
      issubset(tr(5, μ, s), tr(10, μ, s)) for μ in TM_μS, s in 1:S)
    #   - mask not keyed on μ: equal across μ for every sample
    r["masks not all equal across μ"] = all(
      !all(Set(tr(n, TM_μS[1], s)) == Set(tr(n, TM_μS[2], s)) for s in 1:S) for n in (5, 10))
    #   - mask not keyed on the sample: one mask for every sample
    r["masks not all equal across samples"] = all(
      !allequal([Set(tr(n, μ, s)) for s in 1:S]) for n in (5, 10), μ in TM_μS)
  end
  return r
end

"Run `@test` on each named check, printing the names of the failing ones."
function test_checks(r)
  for (name, ok) in sort(collect(r))
    ok || println("FAILED check: ", name)
    @test ok
  end
end

"""
Run `script` with `args`, expecting a refusal before any work (§2 Output labels): a nonzero
exit, an error matching `pattern` (default: one that names the doping counts), and no
output file in `outdir` (the directory itself may exist).
"""
function test_refused(script, outdir, args...; pattern = r"[Dd]oping counts")
  code, out, err = run_script(script, args..., "outdir=$outdir")
  files = isdir(outdir) ? readdir(outdir) : String[]
  ok = code != 0 && occursin(pattern, err) && isempty(files)
  ok || @info "expected a refusal before any output" script args code files err
  @test ok
end

@testset "item 9: count-mode doping in the scripts" begin
  mktempdir() do dir
    @testset "Comptime_prob_tm.jl, four methods, :xy, n = 0, 3, T" begin
      out = joinpath(dir, "ct")
      counts = [0, 3, CT_T]
      code, stdout_, err = run_script("Comptime_prob_tm.jl", CT_ARGS...,
        "magic_counts=$(join(counts, ','))", "outdir=$out")
      @test code == 0
      code == 0 ? test_checks(check_comptime(out; counts, methods = CT_METHODS, mode = :xy)) :
                  @info "Comptime_prob_tm.jl failed" stdout_ err
    end

    @testset "Comptime_prob_tm.jl, CAMPS-PP, :z, n = 3, T" begin
      out = joinpath(dir, "ctz")
      counts = [3, CT_T]
      code, stdout_, err = run_script("Comptime_prob_tm.jl", "N=$CT_N", "samples=$CT_SAMPLES",
        "methods=campspp", "magic_mode=z", "magic_counts=$(join(counts, ','))", "outdir=$out")
      @test code == 0
      code == 0 ? test_checks(check_comptime(out; counts, methods = ["campspp"], mode = :z)) :
                  @info "Comptime_prob_tm.jl (z) failed" stdout_ err
    end

    @testset "Comptime_prob_tm.jl, per-method caps" begin
      out = joinpath(dir, "ctnmax")
      code, stdout_, err = run_script("Comptime_prob_tm.jl", CT_ARGS...,
        "magic_counts=$(join(NMAX_COUNTS, ','))", ("$(m)_nmax=$c" for (m, c) in NMAX_CAPS)...,
        "outdir=$out")
      @test code == 0
      code == 0 ? test_checks(check_nmax(out, stdout_)) :
                  @info "Comptime_prob_tm.jl (caps) failed" stdout_ err
    end

    @testset "Comptime_prob_tm.jl, samples stopped on Pmax" begin
      out = joinpath(dir, "ctpmax")
      code, stdout_, err = run_script("Comptime_prob_tm.jl", PMAX_ARGS..., "outdir=$out")
      @test code == 0
      code == 0 ? test_checks(check_pmax(out, stdout_)) :
                  @info "Comptime_prob_tm.jl (Pmax) failed" stdout_ err
    end

    # A cap below the first count (camps would never run) is refused before any output;
    # campspp's cap (its default) is fine, so the error is camps'
    @testset "refused cap: camps_nmax below every count" begin
      test_refused("Comptime_prob_tm.jl", joinpath(dir, "bad_nmax"), "N=8", "samples=1",
        "methods=campspp,camps", "magic_counts=3,28", "camps_nmax=2";
        pattern = r"camps_nmax = 2 is below every doping count")
    end

    @testset "TransferredMag_varparams.jl" begin
      out = joinpath(dir, "tm")
      code, stdout_, err = run_script("TransferredMag_varparams.jl", "N=$TM_N",
        "Nsamples=$TM_SAMPLES", "μs=$(join(TM_μS, ','))", "magic_counts=$(join(TM_COUNTS, ','))",
        "outdir=$out"; threads = 2)
      @test code == 0
      code == 0 ? test_checks(check_transferredmagnetization(out)) :
                  @info "TransferredMag_varparams.jl failed" stdout_ err
    end

    # §2: 0 ≤ n ≤ T, and the count list strictly increasing; refused before the run starts
    # (review follow-up 8). T = 28 at N = 8 for both scripts (Comptime: N/2 layers; TM: t = N ÷ 2).
    @testset "refused counts: $script $counts" for script in
        ("Comptime_prob_tm.jl", "TransferredMag_varparams.jl"), counts in ("0,29", "0,5,3", "3,3")
      extra = script == "Comptime_prob_tm.jl" ? ("samples=1", "methods=campspp") :
                                                ("Nsamples=1", "μs=0.6")
      test_refused(script, joinpath(dir, "bad_$(script[1:2])_$(replace(counts, ',' => '_'))"),
                   "N=8", "magic_counts=$counts", extra...)
    end
  end
end
