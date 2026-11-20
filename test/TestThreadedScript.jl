# Plan item 2, script level: TransferredMag_varparams.jl, the one script that runs samples in
# `Threads.@threads`, failed at -t 2 with "rotation needs a Hermitian Pauli string" (shared
# QuantumClifford buffers; CampsPP/test/TestThreads.jl tests the cause). Here the plan's
# reproducer runs at -t 1, 2 and 4 in child processes on the root environment:
#   - every run exits 0 and logs the thread count it was given;
#   - the per-sample (`_full.txt`) and averaged (`.txt`) outputs at -t 2 and -t 4 equal those
#     at -t 1 (README: each sample has its own seed, so results don't depend on the thread
#     count). Comment lines must match exactly, numbers to 1e-12. The tolerance dates from when
#     the script set BLAS threads to the Julia thread count, so runs could differ in BLAS
#     rounding (plan, item 2). Since 2026-11-20 (plan, item 10) BLAS/Strided threads follow
#     `blas_threads`, default 1 at every -t, which these runs don't set. The files were
#     byte-identical when this was written (2026-11-19).
# Sensitivity (2026-11-19): with the pre-fix `_apply_on_sites!` patched in memory into a
# process that then includes the script, the reproducer at -t 2 exited 1 in 4 of 4 runs
# ("rotation needs a Hermitian Pauli string").
#
# Slow tier only (STABSANDWICH_SLOW set): three script runs of about 25–30 s each.
#   STABSANDWICH_SLOW=1 julia --project=. -t 2 test/TestThreadedScript.jl   (or test/runtests.jl)

using Test

const ROOT = dirname(@__DIR__)
const SCRIPT = joinpath(ROOT, "TransferredMag_varparams.jl")
const REPRO_ARGS = ["N=8", "Nsamples=3", "μs=0.6", "magic_probs=0,1"]
const SCRIPT_ATOL = 1e-12

"""
Compare two output texts line by line: comment and empty lines exactly, data lines as numbers
to `atol`. Returns `(same, max |difference|)`; `same` is false on any structural difference.
"""
function compare_outputs(a::AbstractString, b::AbstractString; atol = SCRIPT_ATOL)
  la, lb = split(a, '\n'), split(b, '\n')
  length(la) == length(lb) || return (false, Inf)
  maxdiff = 0.0
  for (x, y) in zip(la, lb)
    if isempty(strip(x)) || startswith(x, '#')
      x == y || return (false, Inf)
      continue
    end
    vx, vy = split(x), split(y)
    length(vx) == length(vy) || return (false, Inf)
    nx = tryparse.(Float64, vx); ny = tryparse.(Float64, vy)
    (any(isnothing, nx) || any(isnothing, ny)) && return (false, Inf)
    maxdiff = max(maxdiff, maximum(abs.(nx .- ny)))
  end
  return (maxdiff <= atol, maxdiff)
end

"Data lines (not comments, not empty) of an output text."
ndatalines(text) = count(l -> !isempty(strip(l)) && !startswith(l, '#'), split(text, '\n'))

@testset "compare_outputs (negative checks)" begin
  text = "# header\n\n0\t0.0\n1\t2.0\n"
  @test compare_outputs(text, text) == (true, 0.0)
  @test compare_outputs(text, replace(text, "2.0" => "2.0000000001"))[1] == false   # 1e-10
  @test compare_outputs(text, replace(text, "2.0" => "2.0000000000001"))[1] == true # 1e-13
  @test compare_outputs(text, replace(text, "# header" => "# other"))[1] == false
  @test compare_outputs(text, text * "2\t1.0\n")[1] == false
  @test compare_outputs(text, replace(text, "1\t2.0" => "1\t2.0\t0.0"))[1] == false
end

if get(ENV, "STABSANDWICH_SLOW", "") != ""
@testset "TransferredMag_varparams.jl at -t 1, 2, 4 (reproducer)" begin
  mktempdir() do dir
    outputs = Dict{Int, Dict{String, String}}()
    for n in (1, 2, 4)
      outdir = joinpath(dir, "t$n")
      cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no -t $n $SCRIPT outdir=$outdir $REPRO_ARGS`
      obuf, ebuf = IOBuffer(), IOBuffer()
      proc = run(pipeline(ignorestatus(cmd); stdout = obuf, stderr = ebuf))
      out, err = String(take!(obuf)), String(take!(ebuf))
      @test proc.exitcode == 0
      proc.exitcode == 0 || @info "-t $n failed" out err
      @test occursin("Using $n threads.", out)
      @test !occursin("WARNING", out)
      outputs[n] = Dict(suffix => read(joinpath(outdir, "TMDxy_8_3$suffix.txt"), String)
                        for suffix in ("", "_full") if isfile(joinpath(outdir, "TMDxy_8_3$suffix.txt")))
    end
    # Structure of the -t 1 reference: 2 (p, μ) pairs × t+1 = 5 cycles, 3 samples in _full
    ref = outputs[1]
    @test ndatalines(get(ref, "_full", "")) == 2 * 3 * 5
    @test ndatalines(get(ref, "", "")) == 2 * 5
    for n in (2, 4), suffix in ("", "_full")
      same, maxdiff = compare_outputs(get(outputs[n], suffix, ""), get(ref, suffix, "-"))
      @test same
      same || @info "-t $n vs -t 1, TMDxy_8_3$suffix.txt" maxdiff
    end
  end
end
end
