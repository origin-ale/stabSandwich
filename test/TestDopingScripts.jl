# Post-port item 9, script level: count-mode doping in Comptime_prob_tm.jl and
# TransferredMag_varparams.jl (test/CONVENTIONS.md §2 "Doping": count mode, seeds in the
# scripts, output labels). Each script runs once, small, in a child process on the root
# environment, and its output files are checked against the spec:
#   - Comptime_prob_tm.jl: every sample's doped gates are exactly n XX-YY pairs on n distinct
#     triples, and those triples are the mask of MersenneTwister(hash((rng_seed + i, n)));
#     initial states from seed rng_seed + i; the x column of comptime_prob_tm.txt is n;
#     every block is headed `# n = <integer>, μ = <μ>`.
#   - TransferredMag_varparams.jl: a (μ, sample) has the same initial state for every n
#     (state stream on (μ index, sample) only), masks are independent across n (not nested,
#     as a shared state/doping stream would make them) and across μ, exactly n triples each.
# The spec fixes the Comptime_prob_tm.jl mask seed exactly; for TransferredMag_varparams.jl
# it fixes only what each stream depends on, so those checks are relational.
#
#   julia --project=. -t 2 test/TestDopingScripts.jl   (also run by test/runtests.jl)
# About 70 s, mostly compilation in the two child processes.
#
# STABSANDWICH_SCRIPT_DIR (default: the repository root) points the runs at another copy of
# the scripts; used for the in-memory/scratch mutants of item 9, not by the suite.

using Test
using CampsPP
using CampsPP.DisentangleCAMPS
using Random: MersenneTwister

const ROOT = dirname(@__DIR__)
const SCRIPT_DIR = get(ENV, "STABSANDWICH_SCRIPT_DIR", ROOT)

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
Triples doped in `:xy` mode by the gate list `gs`, or `nothing` unless the gates are exactly
an XX at 3m−2 and a YY at 3m−1 for each listed triple m.
"""
function xy_triples(gs)
  ms = sort(unique([cld(g, 3) for (g, _) in gs]))
  want = sort([[(3m - 2, "XX") for m in ms]; [(3m - 1, "YY") for m in ms]])
  sort([(g, filter(!=('I'), P)) for (g, P) in gs]) == want || return nothing
  return ms
end

# -- Comptime_prob_tm.jl ------------------------------------------------------------------
const CT_N, CT_SAMPLES, CT_SEED, CT_μ = 8, 2, 100, 0.6
const CT_COUNTS = [0, 3]

"Spec checks on a Comptime_prob_tm.jl output directory; returns a Dict of named Bools."
function check_comptime(dir)
  N = CT_N
  _, phases = xxz_circuit(π/4, π/4, N ÷ 2, N)
  T = length(phases) ÷ 3                                  # (N−1)·N/2 = 28
  r = Dict{String, Bool}()

  doped = doped_gates(joinpath(dir, "doped_prob_tm.txt"))
  r["doped file parses"] = doped !== nothing && sort(collect(keys(doped))) == CT_COUNTS &&
                           all(length(doped[n]) == CT_SAMPLES for n in CT_COUNTS)
  if r["doped file parses"]
    ms = Dict((n, s) => xy_triples(doped[n][s]) for n in CT_COUNTS, s in 1:CT_SAMPLES)
    r["2n gates = n XX-YY pairs on n distinct triples"] =
      all(v !== nothing && length(v) == n && all(1 .<= v .<= T) for ((n, _), v) in ms)
    # §2: sample i's mask is doping_mask_count(MersenneTwister(hash((rng_seed + i, n))), …)
    r["mask from hash((rng_seed + i, n))"] = all(
      ms[(n, s)] == findall(doping_mask_count(MersenneTwister(hash((CT_SEED + s, n))), phases, n))
      for n in CT_COUNTS, s in 1:CT_SAMPLES)
  end

  # §2: sample i's initial state from seed rng_seed + i, the same for every n
  states = [split(l) for l in readlines(joinpath(dir, "initstates_prob_tm.txt"))
            if !isempty(l) && !startswith(l, '#')]
  r["initial states from rng_seed + i"] = length(states) == CT_SAMPLES && all(
    parse(Int, f[1]) == s &&
    parse.(Int, f[3:end]) == domainwallstate(MersenneTwister(CT_SEED + s), N, CT_μ)[2]
    for (s, f) in enumerate(states))

  # §2: the x column of comptime_prob_tm.txt (and gctime) is n, written as an integer
  r["x column is n"] = all(begin
      rows = [split(l) for l in readlines(joinpath(dir, f)) if !isempty(l) && !startswith(l, '#')]
      [r[1] for r in rows] == string.(CT_COUNTS)
    end for f in ("comptime_prob_tm.txt", "gctime_prob_tm.txt"))

  # §2: every block label in every output is `# n = <integer>, μ = <μ>`, for each n in order
  files = filter(f -> endswith(f, ".txt") && !endswith(f, "_log.txt"), readdir(dir))
  labels = [(f, l) for f in files for l in label_lines(joinpath(dir, f))]
  r["labels present"] = !isempty(labels)
  r["every label is # n = <int>, μ = 0.6"] = all(occursin(LABEL_RE, l) &&
    match(LABEL_RE, l)[2] == string(CT_μ) for (_, l) in labels)
  r["blocks in n order"] = all(begin
      ns = [b[1] for b in label_blocks(joinpath(dir, f))]
      isempty(ns) || unique(ns) == CT_COUNTS
    end for f in files)
  r["header lists magic_counts"] =
    occursin("magic_counts=0,3", readline(joinpath(dir, "comptime_prob_tm.txt")))
  return r
end

# -- TransferredMag_varparams.jl -----------------------------------------------------------
const TM_N, TM_SAMPLES = 8, 3
const TM_μS = ["0.3", "0.6"]
const TM_COUNTS = [0, 5, 10]

"Spec checks on a TransferredMag_varparams.jl output directory; returns a Dict of named Bools."
function check_transferredmag(dir)
  N = TM_N
  T = (N - 1) * (N ÷ 2)                                    # t = N ÷ 2 layers: 28 triples
  stem = joinpath(dir, "TMDxy_$(N)_$(TM_SAMPLES)")
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

  # Initial states: sample → bitstring, per (n, μ)
  states = Dict((b[1], b[2]) => [split(l) for l in b[3]] for b in label_blocks(stem * "_states.txt"))
  r["states parse"] = all(haskey(states, k) && length(states[k]) == TM_SAMPLES for k in pairs_want) &&
    all(join(q in parse.(Int, f[3:end]) ? '1' : '0' for q in 1:N) == f[2] &&
        parse(Int, f[1]) == s for v in values(states) for (s, f) in enumerate(v))
  if r["states parse"]
    bits(n, μ) = [f[2] for f in states[(n, μ)]]
    # §2: the state stream depends on (μ index, sample) only, so not on n
    r["state identical across n for each (μ, sample)"] =
      all(bits(n, μ) == bits(TM_COUNTS[1], μ) for n in TM_COUNTS, μ in TM_μS)
    r["states differ between samples"] = all(allunique(bits(0, μ)) for μ in TM_μS)
  end

  # Doped triples: sample → set of triples, per (n, μ)
  doped = Dict((b[1], b[2]) => [parse.(Int, split(l)) for l in b[3]]
               for b in label_blocks(stem * "_doped.txt"))
  r["doped parse"] = all(haskey(doped, k) && length(doped[k]) == TM_SAMPLES &&
                         [v[1] for v in doped[k]] == 1:TM_SAMPLES for k in pairs_want)
  if r["doped parse"]
    tr(n, μ, s) = doped[(n, μ)][s][2:end]
    r["exactly n distinct triples in 1:T"] = all(
      length(tr(n, μ, s)) == n && allunique(tr(n, μ, s)) && all(1 .<= tr(n, μ, s) .<= T)
      for n in TM_COUNTS, μ in TM_μS, s in 1:TM_SAMPLES)
    # §2: the mask stream is keyed on (n, μ index, sample). Independent draws at T = 28:
    # equal 5-subsets have probability 1/binomial(28, 5) ≈ 1e-5, nested 5 ⊂ 10 about 2e-3.
    r["masks not nested across n"] = !any(
      issubset(tr(5, μ, s), tr(10, μ, s)) for μ in TM_μS, s in 1:TM_SAMPLES)
    r["masks differ across μ"] = all(
      Set(tr(n, TM_μS[1], s)) != Set(tr(n, TM_μS[2], s)) for n in (5, 10), s in 1:TM_SAMPLES)
    r["masks differ across samples"] = all(
      allunique([Set(tr(n, μ, s)) for s in 1:TM_SAMPLES]) for n in (5, 10), μ in TM_μS)
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

@testset "item 9: count-mode doping in the scripts" begin
  mktempdir() do dir
    @testset "Comptime_prob_tm.jl" begin
      out = joinpath(dir, "ct")
      code, stdout_, err = run_script("Comptime_prob_tm.jl", "N=$CT_N", "magic_counts=0,3",
        "samples=$CT_SAMPLES", "methods=campspp", "outdir=$out")
      @test code == 0
      code == 0 ? test_checks(check_comptime(out)) : @info "Comptime_prob_tm.jl failed" stdout_ err
    end

    @testset "TransferredMag_varparams.jl" begin
      out = joinpath(dir, "tm")
      code, stdout_, err = run_script("TransferredMag_varparams.jl", "N=$TM_N",
        "Nsamples=$TM_SAMPLES", "μs=$(join(TM_μS, ','))", "magic_counts=$(join(TM_COUNTS, ','))",
        "outdir=$out"; threads = 2)
      @test code == 0
      code == 0 ? test_checks(check_transferredmag(out)) : @info "TransferredMag_varparams.jl failed" stdout_ err
    end
  end
end
