# Kicked Ising, script level: Comptime_kim.jl against test/CONVENTIONS.md §2 "Kicked Ising"
# (KIM plan step 2 and its budget revision). The script runs small, in child processes on the
# root environment, and its output files are checked against the spec and the dense reference.
# Runs stop on per-method budgets, checked at layer ends, or on the cap `t_max`; every check
# below follows each run to its own stop layer l (from layers_kim_<m>.txt).
#   1. Budgets that can't bind (`*_Pmax=10^9`, `*_χmax=10^6`): `eps_mode = :uniform`, all four
#      methods, εs = 0, 0.3, sites 1, 4, 6, t_max = 8:
#      - every run stops on the cap: l = t_max, s = 0;
#      - every method's per-run ⟨Z_j⟩ (evs_kim_<m>_<exp>_full.txt) equals the dense ⟨Z_j⟩ of the
#        run's (ε, j), with (i, j) read from runs_kim.txt (what ran, not only what was
#        recorded); run r is (i, j) with r = (i − 1)·|sites| + (position of j);
#      - the profile file's means and errors are those of the runs that recorded each layer;
#        every block of every file is headed `# ε = <ε>` (no μ), in εs order;
#      - comptime_kim.txt, gctime_kim.txt and layers_kim.txt have x column ε and per method
#        the 7 statistics (mean, std. err., q1, median, q3, min, max) over every run of the
#        per-run times (comptime_full_kim_<m>.txt) and stop layers.
#   2. Budgets that bind, chosen from run 1's resource files (PP and CAMPS-PP: numbers of Pauli
#      strings; CAMPS: bond dimension) and, for MPS, from `maxlinkdim` of unstopped
#      `mps_circuit_dynamics` runs on the first l layers (the script writes no MPS bond
#      dimension): same settings otherwise. PP's budget is reached by some run exactly at
#      t_max, which must then report s = 1 (the others bind before the cap).
#      - Each run stops at the first layer whose value is ≥ its method's budget, s = 1, or on
#        the cap with s = 0 if none is (the ε = 0 runs of PP, CAMPS-PP and CAMPS);
#      - the run's resource values are run 1's up to l, so every earlier layer is below the
#        budget; the full files end at l and equal dense up to l;
#      - negative check: some run's value at l equals the budget, so a `>` rule would stop it
#        later.
#   2b. CAMPS-PP, CAMPS and MPS at ε = 0.3 with the cap lowered to the last layer tc at which
#      each of them has a run with a new largest value, and that value as budget: that run
#      reaches its budget exactly on the cap and must report (tc, 1).
#   3. `eps_mode = :random`, samples = 2, all four methods, εs = 0.1, 0.3, sites 5, 2, 4 (not
#      sorted, so a script that sorts them or runs them sorted fails):
#      - the u of each (sample, ε) in eps_kim.txt is
#        2 .* rand(MersenneTwister(hash((rng_seed + i, ε))), N) .- 1, and the u at the two ε
#        differ (independent draws, not one pattern scaled by ε);
#      - CAMPS-PP, CAMPS and MPS (exact, budgets not binding) give the dense ⟨Z_j⟩ with
#        ε_j = ε·u_j;
#      - PP runs with pp_thl = 1e-4 and a binding pp_Pmax, so the samples of a site stop at
#        different layers (at pp_thl = 1e-14 the number of strings does not depend on u): the
#        stop layers are those predicted from `pauliprop_circuit_dynamics` with the same
#        threshold, and the profile rows average only the samples that recorded each layer
#        (NaN mean if none, NaN error below two), up to the deepest layer recorded at that ε.
#      Some ε block mixes runs stopped on the budget and runs on the cap, so layers_kim.txt is
#      checked over both kinds.
#   4. Budgets and thresholds of unselected methods are not refused (`methods=pp` with the
#      other budgets 1 and campspp_thl = 0; `methods=mps` with pp_thl = 0, pp_Pmax = 1,
#      campspp_thl < 0, campspp_Pmax = 1).
#   5. Refused before any output file or directory is written: `:uniform` with samples = 2,
#      `:random` with samples = 0, pp_thl = 0, campspp_thl = 0, εs empty or not strictly
#      increasing (0.3,0 and 0,0), sites empty, outside 1:N (0 and N + 1) or repeated, methods
#      unknown, empty or repeated, an unknown eps_mode, t_max = 0, pp_Pmax = 1,
#      campspp_Pmax = 1, camps_χmax = 1, mps_χmax = 1 (an integer, so the script's own check
#      refuses it, not the integer parsing). A refusal must exit nonzero, name the offending
#      parameter on stderr, and leave `outdir` absent.
#
# Exact regime (§6, 1e-10), runs 1–2: N = 6; campspp_χ = 4, so CAMPS-PP switches to PP
# mid-circuit; mps_thl = 1e-15; pp_thl = campspp_thl = 1e-14 (they must be > 0: PP takes the
# Clifford gates as generic rotations, whose cos(π/2) ≈ 6e-17 branches only a positive
# threshold removes; 1e-14 drops nothing else of consequence at N = 6, the default 1e-10 errs
# by up to 7e-11 here); camps_thl at its default 1e-12.
#
# t_max = 8: from the Néel state ⟨Z_j⟩ = 0 after layers 1 … ⌈N/2⌉ = 3 for every ε
# (CampsPP/test/TestCircuits.jl), so only later layers see ε, u and the sites' circuit. The
# ε = 0 block stays 0 after t = 0 to the cap (it checks the initial state, the cap and the
# bookkeeping). With uniform ε, sites 4 and 6 have equal values (shift by two), so a 4 ↔ 6
# mix-up shows only in the random run, whose sites 2, 4 include a same-parity pair.
#
# Not checked (the spec leaves them open): the threshold tag in the evs_kim_* names (only
# that the full and profile files share it), the PW/PC, non-full and `_peak` resource files,
# campspp_switch_kim.txt and the cross-check between methods.
#
# Output layouts the spec leaves open, as assumed here (KIM plan step 3: "via the existing
# writers"):
#   - the full files evs_kim_<m>_<exp>_full.txt and resources_kim_{NP,bd}_full_<m>.txt as
#     save_full_samples writes them: `# sample <r>` per run (`# run <r>` is accepted too; no
#     warm-up run), then `layer<TAB>value` rows, layers without a value left out;
#   - CAMPS-PP's string counts after the switch under the layer they were recorded at. Where
#     the hand-over entry of the driver's `track` output (the unpropagated observable, 1
#     string) is written, if at all, is not checked: `Comptime_prob_tm.jl` writes it under
#     the last CAMPS layer, and as 1 < every budget it never decides a stop;
#   - runs_kim.txt rows `r i j`; layers_kim_<m>.txt rows `r l s`;
#   - comptime_full_kim_<m>.txt as save_per_sample writes it: rows `r time gctime`;
#   - eps_kim.txt one row per sample, `i u_1 … u_N` (or `u_1 … u_N`), at full (round-trip)
#     precision;
#   - the profile file rows `layer mean_1 err_1 mean_2 err_2 …` in `sites` order;
#   - block labels are recognised as `# <ε|eps|n|p|μ> = …` lines, so other comment lines of
#     the form `# x = …` don't count as labels.
#
# Tiers (KIM plan step 4): the helper and reference checks run in the fast tier; the script
# runs (testset "Comptime_kim.jl (KIM)", about 6 min) only in the slow tier.
#   STABSANDWICH_SLOW=1 julia --project=. -t 2 test/TestKimScript.jl   (or test/runtests.jl)
#
# STABSANDWICH_SCRIPT_DIR (default: the repository root) points the runs at another copy of
# the script; used for the scratch mutants of KIM plan step 4, not by the suite.

using Test
using Random: MersenneTwister
using Statistics: mean, std, quantile

const ROOT = dirname(@__DIR__)
const SCRIPT_DIR = get(ENV, "STABSANDWICH_SCRIPT_DIR", ROOT)
const SCRIPT = "Comptime_kim.jl"

isdefined(@__MODULE__, :DenseReference) ||
  include(joinpath(ROOT, "CampsPP", "test", "DisentangleCAMPS", "DenseReference.jl"))
const DR = DenseReference
include(joinpath(ROOT, "CampsPP", "test", "KimSpec.jl"))

# The drivers, for the references of the stops that the script's files don't give (MPS bond
# dimensions; PP string counts of the random run).
using CampsPP, ITensors, ITensorMPS

"Run `script` in SCRIPT_DIR with `args` in a child process; returns (exit code, stdout, stderr)."
function run_script(script, args...; threads = 1)
  cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no -t $threads $(joinpath(SCRIPT_DIR, script)) $args`
  obuf, ebuf = IOBuffer(), IOBuffer()
  proc = run(pipeline(ignorestatus(cmd); stdout = obuf, stderr = ebuf))
  return proc.exitcode, String(take!(obuf)), String(take!(ebuf))
end

maxabsdiff(a, b) = length(a) == length(b) ? maximum(abs.(a .- b); init = 0.0) : Inf

"Equal up to rounding (values are written at full precision), NaN = NaN."
same(a, b; atol = 1e-13, rtol = 1e-12) = length(a) == length(b) &&
  all(isnan(x) ? isnan(y) : isapprox(x, y; atol, rtol) for (x, y) in zip(a, b))

"The seven summary columns of `xs`: mean, std. err., q1, median, q3, min, max (NaN if empty)."
summary7(xs) = isempty(xs) ? fill(NaN, 7) :
  [mean(xs), std(xs) / sqrt(length(xs)), quantile(xs, [0.25, 0.5, 0.75])..., minimum(xs), maximum(xs)]

# -- Parsing --------------------------------------------------------------------------------

const LABEL_RE = r"^# ε = (\S+?)(, max sample = \d+)?$"

"Lines of `path` that look like a block label (`# ε = …`, or a wrong one: eps, n, p, μ)."
label_lines(path) = filter(l -> occursin(r"^# (ε|eps|n|p|μ) = ", l), readlines(path))

"""
Split `path` into its `# ε = …` blocks: a vector of `(ε string, lines)`, with each block's
other non-empty lines. Lines before the first label are dropped.
"""
function eps_blocks(path)
  blocks = Tuple{String, Vector{String}}[]
  for l in eachline(path)
    m = match(LABEL_RE, l)
    if m !== nothing
      push!(blocks, (String(m[1]), String[]))
    elseif !isempty(blocks) && !isempty(strip(l))
      push!(blocks[end][2], l)
    end
  end
  return blocks
end

"Data rows (comments and blanks dropped) of `lines` as vectors of fields."
datarows(lines) = [split(l) for l in lines if !isempty(strip(l)) && !startswith(l, '#')]

"""
Per-run trajectories of a full file: `[(ε, [(run, layers, values), …]), …]` in file order,
or `nothing` if a `# sample`/`# run` header is out of sequence.
"""
function full_runs(path)
  isfile(path) || return nothing
  out = []
  for (ε, lines) in eps_blocks(path)
    runs = []
    for l in lines
      m = match(r"^# (?:sample|run) (\d+)$", l)
      if m !== nothing
        parse(Int, m[1]) == length(runs) + 1 || return nothing
        push!(runs, (parse(Int, m[1]), Int[], Float64[]))
      elseif !startswith(l, '#') && !isempty(runs)
        f = split(l)
        push!(runs[end][2], parse(Int, f[1])); push!(runs[end][3], parse(Float64, f[2]))
      end
    end
    push!(out, (ε, runs))
  end
  return out
end

"Rows of each `# ε` block of `path` as numbers: `[(ε string, [row, …]), …]`."
block_rows(path) = isfile(path) ?
  [(ε, [parse.(Float64, f) for f in datarows(lines)]) for (ε, lines) in eps_blocks(path)] : []

"The one file in `dir` whose name matches `re` (`nothing` if none or several)."
function onefile(dir, re)
  fs = filter(f -> occursin(re, f), readdir(dir))
  return length(fs) == 1 ? joinpath(dir, only(fs)) : nothing
end

"runs_kim.txt as (r, i, j) triples, each listed once, in file order."
runs_table(dir) = unique([Tuple(parse.(Int, f[1:3]))
                          for f in datarows(readlines(joinpath(dir, "runs_kim.txt")))])

"Spec run table: r = (i − 1)·|sites| + position of j, sample-major."
spec_runs(samples, sites) =
  [((i - 1) * length(sites) + p, i, j) for i in 1:samples for (p, j) in enumerate(sites)]

"§2: u of sample i at ε (`:random`)."
spec_u(i, ε, N; rng_seed = 100) = 2 .* rand(MersenneTwister(hash((rng_seed + i, ε))), N) .- 1

"u per sample from eps_kim.txt: `[(ε string, [u_1, u_2, …]), …]`."
function eps_table(dir, N)
  map(eps_blocks(joinpath(dir, "eps_kim.txt"))) do (ε, lines)
    us = map(enumerate(datarows(lines))) do (i, f)
      length(f) == N + 1 ? (parse(Int, f[1]) == i ? parse.(Float64, f[2:end]) : nothing) :
      length(f) == N ? parse.(Float64, f) : nothing
    end
    (ε, us)
  end
end

"""
Stop of a run whose value at layer l is `vals[l]` (a Dict; layers without a value absent) with
budget B and cap tmax: (first layer ≥ 1 with value ≥ B, 1), or (tmax, 0) if there is none.
"""
function predicted_stop(vals, B, tmax)
  ls = sort([l for (l, v) in vals if 1 <= l <= tmax && v >= B])
  return isempty(ls) ? (tmax, 0) : (first(ls), 1)
end

# -- Checks -----------------------------------------------------------------------------------

const METHODS_ALL = ["campspp", "camps", "mps", "pp"]
const KIM_N, KIM_TMAX = 6, 8
const ATOL = 1e-10
const EXACT = ["mps_thl=1e-15", "pp_thl=1e-14", "campspp_thl=1e-14"]
const HUGE = ["pp_Pmax=1000000000", "campspp_Pmax=1000000000", "camps_χmax=1000000", "mps_χmax=1000000"]
const UNIFORM_εS, UNIFORM_SITES = [0.0, 0.3], [1, 4, 6]
const RANDOM_εS, RANDOM_SITES = [0.1, 0.3], [5, 2, 4]    # not sorted, on purpose
const RANDOM_PP_THL = 1e-4

"""
Checks on output directory `dir` of a run with cap `tmax`: labels, runs table, stop layers
(against `stops(m, ε, r)`, an expected `(l, s)` or `nothing`), per-run ⟨Z_j⟩ against the dense value `dense(i, ε, j)`
(t_max + 1 values) for the methods in `dense_methods`, profile means and errors over the runs
that recorded each layer, summary tables. Returns `(checks, stoplayers)`: a Dict of named
Bools, and `Dict((m, ε) => [(r, l, s), …])` as read.
"""
function check_run(dir; methods, εs, sites, samples, dense, stops, dense_methods = methods,
                   tmax = KIM_TMAX)
  r = Dict{String, Bool}()
  εstr = string.(εs)
  read_stops = Dict{Tuple{String, Float64}, Vector{NTuple{3, Int}}}()

  # Labels: every label line of every file is `# ε = <ε>` (no μ), blocks in εs order
  files = filter(f -> endswith(f, ".txt") && !endswith(f, "_log.txt"), readdir(dir))
  labels = [(f, l) for f in files for l in label_lines(joinpath(dir, f))]
  r["some files have ε blocks"] = length(unique(first.(labels))) >= 3 * length(methods)
  r["every label is # ε = <ε> (no μ)"] = all(occursin(LABEL_RE, l) for (_, l) in labels)
  r["no label mentions μ"] = !any(occursin("μ", l) for (_, l) in labels)
  r["blocks in εs order, values parse to εs"] = all(begin
      bs = [b[1] for b in eps_blocks(joinpath(dir, f))]
      isempty(bs) || (bs == εstr && parse.(Float64, bs) == εs)
    end for f in files)

  # Runs: sample-major numbering
  runs = runs_table(dir)
  R = samples * length(sites)
  r["runs_kim.txt is sample-major"] = runs == spec_runs(samples, sites)
  site_of = Dict(rr => j for (rr, _, j) in runs)
  sample_of = Dict(rr => i for (rr, i, _) in runs)

  ls_of = Dict{Tuple{String, Float64}, Vector{Int}}()     # stop layers, per (method, ε)
  for m in methods
    # Stop layers: one `r l s` row per run, l ≤ t_max, s ∈ {0, 1}, s = 0 only on the cap
    lb = block_rows(joinpath(dir, "layers_kim_$m.txt"))
    r["$m: layers_kim_$m.txt, one block per ε, rows r l s"] = [b[1] for b in lb] == εstr &&
      all([Int(row[1]) for row in rows] == 1:R && all(length(row) == 3 for row in rows) for (_, rows) in lb)
    r["$m: layers_kim_$m.txt, one block per ε, rows r l s"] || continue
    for ((_, rows), ε) in zip(lb, εs)
      read_stops[(m, ε)] = [Tuple(Int.(row)) for row in rows]
      ls_of[(m, ε)] = [Int(row[2]) for row in rows]
    end
    r["$m: 1 ≤ l ≤ t_max, s ∈ {0, 1}, s = 0 only at l = t_max"] = all(
      1 <= l <= tmax && s in (0, 1) && (s == 1 || l == tmax) for v in values(read_stops) for (_, l, s) in v)
    r["$m: stop layers as predicted"] = all(begin
        want = stops(m, ε, rr)
        want === nothing || want == (l, s)
      end for ε in εs for (rr, l, s) in read_stops[(m, ε)])

    full = onefile(dir, Regex("^evs_kim_$(m)_[^_]+_full\\.txt\$"))
    prof = onefile(dir, Regex("^evs_kim_$(m)_[^_]+\\.txt\$"))
    r["$m: one full file, one profile file"] = full !== nothing && prof !== nothing
    (full === nothing || prof === nothing) && continue
    # The threshold tag (left open by the spec) is the same in both names
    r["$m: full and profile files have the same tag"] =
      basename(full) == replace(basename(prof), r"\.txt$" => "_full.txt")
    fr = full_runs(full)
    r["$m: full file parses, one block per ε, every run"] = fr !== nothing &&
      [b[1] for b in fr] == εstr && all(length(b[2]) == R for b in fr)
    r["$m: full file parses, one block per ε, every run"] || continue

    # Each run is recorded at layers 0 … l, its stop layer
    r["$m: each run recorded at layers 0 … l"] = all(
      layers == 0:read_stops[(m, ε)][k][2] for ((_, rs), ε) in zip(fr, εs) for (k, (_, layers, _)) in enumerate(rs))

    # What ran: each run's ⟨Z_j⟩ = dense for its (i, j) from runs_kim.txt, up to its stop
    if m in dense_methods
      errs = [maxabsdiff(vals, dense(sample_of[rr], ε, site_of[rr])[1:length(vals)])
              for ((_, rs), ε) in zip(fr, εs) for (rr, _, vals) in rs]
      ok = !isempty(errs) && maximum(errs) <= ATOL
      ok || println("  $m: max |⟨Z_j⟩ − dense| = ", maximum(errs; init = -1.0))
      r["$m: per-run ⟨Z_j⟩ = dense for the run's (sample, ε, site)"] = ok
    end

    # Profile: per layer up to the deepest one recorded at this ε, mean and std. err. over the
    # samples of each site that recorded the layer (NaN mean if none, NaN error below two)
    pb = block_rows(prof)
    r["$m: profile one block per ε"] = [b[1] for b in pb] == εstr
    r["$m: profile one block per ε"] || continue
    r["$m: profile = mean, std. err. over the runs that reached each layer"] = all(begin
        deepest = maximum(length(vals) - 1 for (_, _, vals) in rs)
        length(rows) == deepest + 1 && all(begin
            row = rows[l + 1]
            length(row) == 1 + 2 * length(sites) && row[1] == l && all(begin
                xs = [vals[l + 1] for (rr, _, vals) in rs if site_of[rr] == j && length(vals) > l]
                same([row[2p], row[2p + 1]],
                     [isempty(xs) ? NaN : mean(xs), length(xs) < 2 ? NaN : std(xs) / sqrt(length(xs))])
              end for (p, j) in enumerate(sites))
          end for l in 0:deepest)
      end for ((_, rows), (_, rs)) in zip(pb, fr))
  end

  # Summary tables: x column ε, then per method the 7 statistics over every run, of the
  # per-run times (comptime, gctime) and of the stop layers (layers_kim.txt)
  times = Dict(m => Dict(parse(Float64, ε) => rows for (ε, rows) in
                         block_rows(joinpath(dir, "comptime_full_kim_$m.txt"))) for m in methods)
  r["comptime_full_kim_<m>.txt: every run, r time gctime"] = all(
    haskey(times[m], ε) && [Int(row[1]) for row in times[m][ε]] == 1:R &&
    all(length(row) == 3 for row in times[m][ε]) for m in methods, ε in εs)
  for (f, stat) in (("comptime_kim.txt", (m, ε) -> [row[2] for row in times[m][ε]]),
                    ("gctime_kim.txt", (m, ε) -> [row[3] for row in times[m][ε]]),
                    ("layers_kim.txt", (m, ε) -> Float64.(ls_of[(m, ε)])))
    rows = isfile(joinpath(dir, f)) ? [parse.(Float64, row) for row in datarows(readlines(joinpath(dir, f)))] : []
    shape = [row[1] for row in rows] == εs && all(length(row) == 1 + 7 * length(methods) for row in rows)
    r["$f: x column is ε, 1 + 7·$(length(methods)) fields"] = shape
    shape && r["comptime_full_kim_<m>.txt: every run, r time gctime"] && all(haskey(ls_of, (m, ε)) for m in methods, ε in εs) || continue
    r["$f: per method, the 7 statistics over every run"] = all(
      same(row[1 + 7 * (k - 1) .+ (1:7)], summary7(stat(m, ε)))
      for (row, ε) in zip(rows, εs) for (k, m) in enumerate(methods))
  end
  return r, read_stops
end

"Run `@test` on each named check, printing the names of the failing ones."
function test_checks(r)
  for (name, ok) in sort(collect(r))
    ok || println("FAILED check: ", name)
    @test ok
  end
end

"""
Run the script with `args`, expecting a refusal before any output: a nonzero exit, stderr
matching `pattern` (it names the offending parameter), and `outdir` not created.
"""
function test_refused(outdir, pattern, args...)
  code, out, err = run_script(SCRIPT, args..., "outdir=$outdir")
  ok = code != 0 && occursin(pattern, err) && !ispath(outdir)
  ok || @info "expected a refusal before any output" args code ispath(outdir) err
  @test ok
end

# -- Self-checks of the helpers (no script needed) ---------------------------------------------

@testset "TestKimScript helpers" begin
  @test spec_runs(2, [2, 5]) == [(1, 1, 2), (2, 1, 5), (3, 2, 2), (4, 2, 5)]
  @test spec_runs(1, [1, 4, 6]) == [(1, 1, 1), (2, 1, 4), (3, 1, 6)]
  # §2: u ∈ [−1, 1]^N, reproducible, different per sample and per ε
  u = spec_u(1, 0.1, 6)
  @test length(u) == 6 && all(-1 .<= u .<= 1) && u == spec_u(1, 0.1, 6)
  @test u != spec_u(2, 0.1, 6) && u != spec_u(1, 0.3, 6)
  # Stops: first layer with value ≥ B; the cap with s = 0 if none; layer 0 never stops
  v = Dict(0 => 5, 1 => 2, 2 => 4, 3 => 4, 4 => 9)
  @test predicted_stop(v, 4, 8) == (2, 1)
  @test predicted_stop(v, 5, 8) == (4, 1)
  @test predicted_stop(v, 10, 8) == (8, 0)
  @test predicted_stop(v, 9, 4) == (4, 1)          # reached on the cap: s = 1
  @test predicted_stop(Dict(3 => 7), 2, 8) == (3, 1)   # values from layer 3 only (CAMPS-PP)
  @test summary7([1.0, 2.0, 3.0, 4.0])[[1, 3, 4, 5, 6, 7]] == [2.5, 1.75, 2.5, 3.25, 1.0, 4.0]
  @test isnan(summary7([2.0])[2]) && all(isnan, summary7(Float64[]))
  mktempdir() do dir
    path = joinpath(dir, "f.txt")
    write(path, "# N=6  \n# Cycle x\n\n\n# ε = 0.0\n# sample 1\n0\t1.0\n1\t0.5\n\n# sample 2\n0\t-1.0\n1\t0.25\n\n\n# ε = 0.3, max sample = 2\n# sample 1\n0\t2.0\n\n\n")
    @test [b[1] for b in eps_blocks(path)] == ["0.0", "0.3"]
    fr = full_runs(path)
    @test fr[1][2] == [(1, [0, 1], [1.0, 0.5]), (2, [0, 1], [-1.0, 0.25])]
    @test fr[2][2] == [(1, [0], [2.0])]
    @test length(label_lines(path)) == 2
    write(path, "# r = run\n# ε = 0.1\n1 4 1\n2 8 0\n\n\n")
    @test length(label_lines(path)) == 1           # `# r = …` is not a block label
    @test block_rows(path) == [("0.1", [[1.0, 4.0, 1.0], [2.0, 8.0, 0.0]])]
  end
end

@testset "references of the script runs (no script needed)" begin
  N, t = KIM_N, KIM_TMAX
  # Uniform runs: the ε = 0.3 profile is not 0 after t = 0, and differs between neighbouring
  # sites
  p = kim_dense_profile(0.3, t, N)
  @test maximum(abs, p[2:end, :]) > 0.05
  @test maximum(abs.(p[2:end, 1] .- p[2:end, 2])) > 0.05
  # Random run: another sample's u, or the other ε's u (a u shared across ε), gives other
  # values at the sites run; so do sites 2 and 4 (a same-parity swap) and sites 5 and 2 (a
  # script that runs the sites sorted puts 2 first)
  prof(i, ε; uε = ε) = kim_dense_profile(ε .* spec_u(i, uε, N), t, N)[:, RANDOM_SITES]
  col(j) = findfirst(==(j), RANDOM_SITES)
  @test !issorted(RANDOM_SITES)
  @test all(maxabsdiff(prof(1, ε), prof(2, ε)) > 1e-3 for ε in RANDOM_εS)
  @test all(maxabsdiff(prof(i, 0.3), prof(i, 0.3; uε = 0.1)) > 1e-3 for i in 1:2)
  @test all(maxabsdiff(prof(i, ε)[:, col(2)], prof(i, ε)[:, col(4)]) > 1e-3 for i in 1:2, ε in RANDOM_εS)
  @test all(maxabsdiff(prof(i, ε)[:, col(5)], prof(i, ε)[:, col(2)]) > 1e-3 for i in 1:2, ε in RANDOM_εS)
end

if get(ENV, "STABSANDWICH_SLOW", "") != ""
@testset "Comptime_kim.jl (KIM)" begin
  N, tmax = KIM_N, KIM_TMAX
  L = kim_spec_L(N)
  le = [L * l for l in 1:tmax]
  neel = kim_neel(N)
  # MPS bond dimension after layers 1 … t_max (unstopped runs on the first l layers)
  function mps_bonddims(εv)
    g, ϕ = kicked_ising_circuit(εv, tmax, N)
    Dict(l => maxlinkdim(first(mps_circuit_dynamics(
                MPS(ComplexF64, siteinds("Qubit", N), [j in neel ? "1" : "0" for j in 1:N]),
                g[1:L*l], ϕ[1:L*l], 1e-15, local_z(N, 1), nothing; layer_ends = [L * l])))
         for l in 1:tmax)
  end
  # PP string counts after layers 0 … t_max, at `thl`, of the run (εv, site j)
  function pp_counts(εv, j, thl)
    g, ϕ = kicked_ising_circuit(εv, tmax, N)
    _, evs, nt = pauliprop_circuit_dynamics(neel, g, ϕ, thl, typemax(Int), local_z(N, j);
      layer_ends = le, track = true)
    @assert length(nt) == length(evs) == tmax + 1
    Dict(l => nt[l + 1] for l in 0:tmax)
  end

  mktempdir() do dir
    uniform_args = ["N=$N", "t_max=$tmax", "methods=$(join(METHODS_ALL, ','))",
      "εs=$(join(UNIFORM_εS, ','))", "sites=$(join(UNIFORM_SITES, ','))", "campspp_χ=4", EXACT...]
    uprofiles = Dict(ε => kim_dense_profile(ε, tmax, N) for ε in UNIFORM_εS)
    udense = (i, ε, j) -> uprofiles[ε][:, j]
    out1 = joinpath(dir, "uniform")

    @testset "1. uniform ε, budgets not binding" begin
      code, stdout_, err = run_script(SCRIPT, uniform_args..., HUGE..., "outdir=$out1")
      @test code == 0
      if code == 0
        checks, _ = check_run(out1; methods = METHODS_ALL, εs = UNIFORM_εS, sites = UNIFORM_SITES,
          samples = 1, dense = udense, stops = (m, ε, r) -> (tmax, 0))
        test_checks(checks)
        @test !isfile(joinpath(out1, "eps_kim.txt"))       # :random only
      else
        @info "Comptime_kim.jl (uniform) failed" stdout_ err
      end
    end

    @testset "2. uniform ε, budgets that bind" begin
      # References per (method, ε, run): value at each layer, from run 1's resource files
      # (PP, CAMPS-PP: strings; CAMPS: bond dimension) and from the MPS driver
      res_file = Dict("pp" => "resources_kim_NP_full_pp.txt",
                      "campspp" => "resources_kim_NP_full_campspp.txt",
                      "camps" => "resources_kim_bd_full_camps.txt")
      ref = Dict{Tuple{String, Float64, Int}, Dict{Int, Int}}()
      R = length(UNIFORM_SITES)
      ok_files = isdir(out1) && all(begin
          fr = full_runs(joinpath(out1, f))
          fr !== nothing && [b[1] for b in fr] == string.(UNIFORM_εS) && all(length(b[2]) == R for b in fr)
        end for f in values(res_file))
      @test ok_files
      if ok_files
        for (m, f) in res_file, ((_, rs), ε) in zip(full_runs(joinpath(out1, f)), UNIFORM_εS)
          for (rr, layers, vals) in rs
            ref[(m, ε, rr)] = Dict(l => round(Int, v) for (l, v) in zip(layers, vals))
          end
        end
        for ε in UNIFORM_εS
          bd = mps_bonddims(fill(ε, N))
          for rr in 1:R
            ref[("mps", ε, rr)] = bd
          end
        end
        # Budget per method: the smallest over the ε = 0.3 runs of the largest value before
        # the cap, so every ε = 0.3 run reaches it before t_max. PP instead gets the smallest
        # value at layer t_max, so some run reaches its budget exactly on the cap, where the
        # spec wants s = 1 (a script that sets s = (l < t_max) fails there).
        budget = Dict(m => minimum(maximum(v for (l, v) in ref[(m, 0.3, rr)] if 1 <= l <= tmax - 1)
                                   for rr in 1:R) for m in METHODS_ALL)
        budget["pp"] = minimum(ref[("pp", 0.3, rr)][tmax] for rr in 1:R)
        want(m, ε, rr) = predicted_stop(ref[(m, ε, rr)], budget[m], tmax)
        @testset "premises" begin
          @test all(>=(2), values(budget))                       # the script refuses < 2
          @test all(want(m, 0.3, rr)[1] < tmax && want(m, 0.3, rr)[2] == 1
                    for m in ("campspp", "camps", "mps"), rr in 1:R)
          @test all(want("pp", 0.3, rr)[2] == 1 for rr in 1:R)
          @test any(want("pp", 0.3, rr) == (tmax, 1) for rr in 1:R)    # budget reached on the cap
          # ε = 0: PP keeps one string, CAMPS-PP never switches, CAMPS stays at bond dim. 1
          @test all(want(m, 0.0, rr) == (tmax, 0) for m in ("pp", "campspp", "camps"), rr in 1:R)
          # Negative: a value equal to the budget at the stop layer (a `>` rule stops later)
          @test all(any(ref[(m, 0.3, rr)][want(m, 0.3, rr)[1]] == budget[m] for rr in 1:R) for m in METHODS_ALL)
          # … and the budget does stop some run after its first layer with a value
          @test all(any(want(m, 0.3, rr)[1] > minimum(l for l in keys(ref[(m, 0.3, rr)]) if l >= 1)
                        for rr in 1:R) for m in METHODS_ALL)
        end

        out2 = joinpath(dir, "binding")
        code, stdout_, err = run_script(SCRIPT, uniform_args...,
          "pp_Pmax=$(budget["pp"])", "campspp_Pmax=$(budget["campspp"])",
          "camps_χmax=$(budget["camps"])", "mps_χmax=$(budget["mps"])", "outdir=$out2")
        @test code == 0
        if code == 0
          checks, _ = check_run(out2; methods = METHODS_ALL, εs = UNIFORM_εS, sites = UNIFORM_SITES,
            samples = 1, dense = udense, stops = want)
          # The resource values are run 1's up to the stop layer, so each earlier layer is
          # below the budget
          for (m, f) in res_file
            fr = full_runs(joinpath(out2, f))
            checks["$m: resource values = unstopped run's, up to the stop layer"] = fr !== nothing &&
              [b[1] for b in fr] == string.(UNIFORM_εS) && all(begin
                l = want(m, ε, rr)[1]
                Dict(k => round(Int, v) for (k, v) in zip(layers, vals)) ==
                  filter(kv -> kv[1] <= l, ref[(m, ε, rr)])
              end for ((_, rs), ε) in zip(fr, UNIFORM_εS) for (rr, layers, vals) in rs)
          end
          test_checks(checks)
        else
          @info "Comptime_kim.jl (binding budgets) failed" stdout_ err
        end

        # 2b. CAMPS-PP, CAMPS and MPS reach their budgets exactly on the cap, so s = 1 at
        # l = t_max for them too (run 2 checks that for PP only). The cap tc is the last layer
        # at which each of the three has some ε = 0.3 run with a new largest value (≥ 2); that
        # value is its budget. Same settings as run 1 up to tc, so run 1's values apply.
        @testset "2b. budgets reached on the cap" begin
          cap_ms = ["campspp", "camps", "mps"]
          newmax(v, l) = haskey(v, l) && v[l] >= 2 && all(v[k] < v[l] for k in keys(v) if 1 <= k < l)
          tcs = [l for l in 2:tmax if all(any(newmax(ref[(m, 0.3, rr)], l) for rr in 1:R) for m in cap_ms)]
          @test !isempty(tcs)                                      # premise
          if !isempty(tcs)
            tc = maximum(tcs)
            bcap = Dict(m => minimum(ref[(m, 0.3, rr)][tc] for rr in 1:R if newmax(ref[(m, 0.3, rr)], tc))
                        for m in cap_ms)
            wantc(m, ε, rr) = predicted_stop(filter(kv -> kv[1] <= tc, ref[(m, ε, rr)]), bcap[m], tc)
            # Premise: for each method some run reaches its budget at tc (a script with
            # s = (l < t_max) reports s = 0 there)
            @test all(any(wantc(m, 0.3, rr) == (tc, 1) for rr in 1:R) for m in cap_ms)
            out2b = joinpath(dir, "on_cap")
            code, stdout_, err = run_script(SCRIPT, "N=$N", "t_max=$tc", "methods=$(join(cap_ms, ','))",
              "εs=0.3", "sites=$(join(UNIFORM_SITES, ','))", "campspp_χ=4", EXACT...,
              "campspp_Pmax=$(bcap["campspp"])", "camps_χmax=$(bcap["camps"])",
              "mps_χmax=$(bcap["mps"])", "outdir=$out2b")
            @test code == 0
            if code == 0
              checks, _ = check_run(out2b; methods = cap_ms, εs = [0.3], sites = UNIFORM_SITES,
                samples = 1, dense = udense, stops = wantc, tmax = tc)
              test_checks(checks)
            else
              @info "Comptime_kim.jl (budgets on the cap) failed" stdout_ err
            end
          end
        end
      end
    end

    @testset "3. random ε, samples = 2, ragged profile" begin
      εs, sites, samples = RANDOM_εS, RANDOM_SITES, 2
      rprofiles = Dict((i, ε) => kim_dense_profile(ε .* spec_u(i, ε, N), tmax, N)
                       for i in 1:samples, ε in εs)
      # PP's stops at pp_thl = 1e-4: string counts from the driver, per (ε, sample, site)
      counts = Dict((ε, i, j) => pp_counts(ε .* spec_u(i, ε, N), j, RANDOM_PP_THL)
                    for ε in εs, i in 1:samples, j in sites)
      stop_of(B, ε, i, j) = predicted_stop(counts[(ε, i, j)], B, tmax)
      ragged(B) = count(stop_of(B, ε, 1, j)[1] != stop_of(B, ε, 2, j)[1] for ε in εs, j in sites)
      mixed(B) = any(length(unique(stop_of(B, ε, i, j)[2] for i in 1:samples, j in sites)) == 2 for ε in εs)
      cands = sort(unique(v for c in values(counts) for (l, v) in c if 1 <= l <= tmax && v >= 2))
      B = argmax(b -> (mixed(b), ragged(b), -b), cands)
      @testset "premises" begin
        @test ragged(B) >= 2              # the samples of a site stop at different layers
        @test mixed(B)                    # an ε block has budget stops and cap stops
      end
      want(m, ε, rr) = m != "pp" ? (tmax, 0) :
        stop_of(B, ε, cld(rr, length(sites)), sites[mod1(rr, length(sites))])

      out3 = joinpath(dir, "random")
      code, stdout_, err = run_script(SCRIPT, "N=$N", "t_max=$tmax", "methods=$(join(METHODS_ALL, ','))",
        "eps_mode=random", "samples=$samples", "εs=$(join(εs, ','))", "sites=$(join(sites, ','))",
        "campspp_χ=4", "campspp_thl=1e-14", "mps_thl=1e-15", "pp_thl=$RANDOM_PP_THL", "pp_Pmax=$B",
        "campspp_Pmax=1000000000", "camps_χmax=1000000", "mps_χmax=1000000", "outdir=$out3")
      @test code == 0
      if code == 0
        checks, _ = check_run(out3; methods = METHODS_ALL, εs, sites, samples,
          dense = (i, ε, j) -> rprofiles[(i, ε)][:, j], stops = want,
          dense_methods = ["campspp", "camps", "mps"])
        table = isfile(joinpath(out3, "eps_kim.txt")) ? eps_table(out3, N) : []
        checks["eps_kim.txt: one block per ε, u for each sample"] =
          [b[1] for b in table] == string.(εs) &&
          all(length(us) == samples && all(!isnothing, us) for (_, us) in table)
        if checks["eps_kim.txt: one block per ε, u for each sample"]
          checks["eps_kim.txt: u = 2 .* rand(MersenneTwister(hash((rng_seed + i, ε))), N) .- 1"] =
            all(maxabsdiff(us[i], spec_u(i, ε, N)) <= 1e-15 for ((_, us), ε) in zip(table, εs)
                for i in 1:samples)
          checks["eps_kim.txt: u differ across ε (independent draws)"] =
            all(maxabsdiff(table[1][2][i], table[2][2][i]) > 0.1 for i in 1:samples)
          checks["eps_kim.txt: u differ across samples"] =
            all(maxabsdiff(us[1], us[2]) > 0.1 for (_, us) in table)
        end
        test_checks(checks)
      else
        @info "Comptime_kim.jl (random) failed" stdout_ err
      end
    end

    @testset "4. settings of an unselected method are not refused: $m" for (k, (m, args)) in enumerate([
        ("pp", ["pp_thl=1e-14", "campspp_thl=0", "campspp_Pmax=1", "camps_χmax=1", "mps_χmax=1"]),
        ("mps", ["mps_thl=1e-15", "pp_thl=0", "pp_Pmax=1", "campspp_thl=-1", "campspp_Pmax=1"])])
      out4 = joinpath(dir, "unselected_$k")
      code, stdout_, err = run_script(SCRIPT, "N=$N", "t_max=1", "methods=$m", "εs=0.1",
        "sites=1", args..., "outdir=$out4")
      ok = code == 0 && isfile(joinpath(out4, "layers_kim_$m.txt"))
      ok || @info "Comptime_kim.jl (unselected settings) failed" stdout_ err
      @test ok
    end

    @testset "5. refused before any output: $name" for (k, (name, pattern, args)) in enumerate([
        ("uniform with samples = 2", r"samples", ["eps_mode=uniform", "samples=2"]),
        ("pp_thl = 0", r"(?<!campspp_)pp_thl", ["methods=pp", "pp_thl=0"]),
        ("campspp_thl = 0", r"campspp_thl", ["methods=campspp", "campspp_thl=0"]),
        ("εs = 0.3, 0", r"ε", ["εs=0.3,0"]),
        ("εs = 0, 0", r"ε", ["εs=0,0"]),
        ("site 0", r"site", ["sites=0,2"]),
        ("site N + 1", r"site", ["sites=2,$(N + 1)"]),
        ("repeated site", r"site", ["sites=2,2"]),
        ("unknown method", r"[Mm]ethod", ["methods=campspp,foo"]),
        ("empty methods", r"[Mm]ethod", ["methods="]),
        ("repeated method", r"[Mm]ethod", ["methods=pp,pp"]),
        ("samples = 0 (random)", r"samples", ["eps_mode=random", "samples=0"]),
        ("empty εs", r"ε", ["εs="]),
        ("empty sites", r"site", ["sites="]),
        ("unknown eps_mode", r"eps_mode", ["eps_mode=gaussian"]),
        ("t_max = 0", r"t_max", ["t_max=0"]),
        ("pp_Pmax = 1", r"(?<!campspp_)pp_Pmax", ["methods=pp", "pp_Pmax=1"]),
        ("campspp_Pmax = 1", r"campspp_Pmax", ["methods=campspp", "campspp_Pmax=1"]),
        ("camps_χmax = 1", r"camps_χmax", ["methods=camps", "camps_χmax=1"]),
        ("mps_χmax = 1", r"mps_χmax", ["methods=mps", "mps_χmax=1"])])
      # A neutral directory name, so a script that echoes its outdir cannot match `pattern`
      test_refused(joinpath(dir, "refused_$k"), pattern, "N=$N",
                   (any(startswith("t_max="), args) ? () : ("t_max=2",))..., args...)
    end
  end
end
end
