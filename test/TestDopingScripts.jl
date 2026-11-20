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
# The spec fixes the Comptime_prob_tm.jl seeds exactly; for TransferredMag_varparams.jl it
# fixes only what each stream depends on, so those checks are relational. Each relational
# check asserts that the signature of the wrong seeding (e.g. "nested for every sample")
# does not hold for every sample; the comment at each check bounds the probability that
# independent draws show the signature everywhere, so none of them passes by luck.
#
#   julia --project=. -t 2 test/TestDopingScripts.jl   (also run by test/runtests.jl)
# About 3 min: three script runs of 30–60 s (mostly compilation) and six refused runs of 5 s.
#
# STABSANDWICH_SCRIPT_DIR (default: the repository root) points the runs at another copy of
# the scripts; used for the scratch mutants of item 9, not by the suite.

using Test
using CampsPP
using CampsPP.DisentangleCAMPS
using Random: MersenneTwister

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

  # §2: sample i's mask is doping_mask_count(MersenneTwister(hash((rng_seed + i, n))), …)
  _, phases = xxz_circuit(π/4, π/4, t, N)
  spec_mask(n, s) = findall(doping_mask_count(MersenneTwister(hash((CT_SEED + s, n))), phases, n))
  spec_ones(s) = domainwallstate(MersenneTwister(CT_SEED + s), N, CT_μ)[2]

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
  dense = Dict((n, s) => dense_tm(N, t, spec_ones(s),
                 spec_phases(N, t, spec_mask(n, s), mode; ϕ_xy = π/4, ϕ_z = π/4,
                             magicphase = CT_MAGIC[mode]))
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
exit, an error that names the doping counts, and no output file in `outdir` (the directory
itself may exist).
"""
function test_refused(script, outdir, args...)
  code, out, err = run_script(script, args..., "outdir=$outdir")
  files = isdir(outdir) ? readdir(outdir) : String[]
  ok = code != 0 && occursin(r"[Dd]oping counts", err) && isempty(files)
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
