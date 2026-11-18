# Phase 7: plain-text progress output (`CampsPP/src/Progress.jl`) and its use in the
# evolution loops.
#
# Spec (plan, "Plain-text progress output", and the Phase 7 decisions): one complete,
# flushed line per report, `[yyyy-mm-dd HH:MM:SS] <indent>desc i/n (pct%) key=v …`, with
# `key => (v, vmax)` shown as `key=v (max vmax)`; no `\r` or ANSI codes. A line is due at
# the last step, when the fraction done reaches the next multiple of `frac` (default 10 %),
# or `dt` seconds (default 30) after the previous line; `frac = 0` prints every step. The
# percentage is floored. `log_done!` prints the final line unless that step was printed.
#
# The milestone rule is checked against an integer oracle (`oracle_steps`), written from the
# spec, not from the Float64 code: with frac = a/b, the number of milestones reached at step
# i of n is fld(i·b, n·a). The `dt` trigger is switched off (`dt = Inf`) wherever the
# printed steps are compared exactly, and tested on its own with `dt = 0` and with a sleep.
#
# The loop tests at the end reuse `doped_xxz`, `spec_layer_ends`, `test_observable` and
# `basis_mps` from TestTimeEvolution.jl (same suite module, included before this file).

using Dates

const TS_FORMAT = dateformat"yyyy-mm-dd HH:MM:SS"

"Split a printed line into its timestamp and the text after `] `; `nothing` if malformed."
function split_line(line)
  m = match(r"^\[(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)\] (.*)$"s, line)
  isnothing(m) && return nothing
  return DateTime(m[1], TS_FORMAT), String(m[2])
end

"Text after the timestamp of each line printed to `io` so far (and empties `io`)."
function bodies(io::IOBuffer)
  out = String(take!(io))
  isempty(out) && return String[]
  @assert endswith(out, '\n')
  return [split_line(l)[2] for l in split(chop(out), '\n')]
end

"Step printed in a body `desc i/n (pct%) …`."
printed_step(body) = parse(Int, match(r" (\d+)/\d+ \(\d+%\)", body)[1])

"""
Steps at which a line is due when `log_progress!` is called at each of `calls` (increasing)
for a loop of `n` steps, with `frac` given exactly as a Rational and no `dt` trigger.
"""
function oracle_steps(n, frac::Rational, calls)
  printed = Int[]
  last = -1; reached = 0
  milestones(i) = fld(i * denominator(frac), n * numerator(frac))
  for i in calls
    i == last && continue
    if i >= n || frac == 0 || milestones(i) > reached
      push!(printed, i)
      last = i
      (frac > 0 && n > 0) && (reached = max(reached, milestones(i)))
    end
  end
  return printed
end

"Steps printed by a `ProgressLog` when `log_progress!` is called at each of `calls`."
function logged_steps(n, frac, calls)
  io = IOBuffer()
  p = ProgressLog("x", n; io, frac = Float64(frac), dt = Inf)   # frac: Rational or Float64
  returned = [i for i in calls if log_progress!(p, i)]
  steps = printed_step.(bodies(io))
  @assert steps == returned          # log_progress! returns whether it printed
  return steps
end

"Run `f()` with stdout redirected to a file; return its result and the printed lines."
function capture_stdout(f)
  mktemp() do path, io
    result = redirect_stdout(f, io)
    flush(io)
    return result, readlines(path)
  end
end

"""
Check the lines one loop printed: format `  desc i/n (pct%) keys…` with the given value
keys, n = `n`, the floored percentage, and the steps the oracle expects for calls at
`calls` followed by the final line at `last(calls)` (or `final` if the loop ran no step).
"""
function check_loop_lines(lines, desc, n, keys, calls; final = last(calls))
  isempty(lines) && return false
  steps = Int[]
  for l in lines
    sl = split_line(l)
    isnothing(sl) && return false
    m = match(Regex("^  \\Q$desc\\E (\\d+)/$n \\((\\d+)%\\)(.*)\$"), sl[2])
    isnothing(m) && return false
    i = parse(Int, m[1])
    parse(Int, m[2]) == (n == 0 ? 100 : fld(100 * i, n)) || return false
    # keys in order, after dropping the " (max …)" of tuple values
    vals = replace(m[3], r" \(max [^)]*\)" => "")
    [k[1] for k in eachmatch(r" ([^=]+?)=\S*", vals)] == keys || return false
    push!(steps, i)
  end
  want = oracle_steps(n, 1//10, calls)
  (isempty(want) || last(want) != final) && push!(want, final)
  return steps == want
end

@testset "Progress.jl" begin
  @testset "line format" begin
    io = IOBuffer()
    t0 = floor(now(), Second)
    p = ProgressLog("CAMPS gate", 10; io, frac = 0, dt = Inf)
    @test log_progress!(p, 3, "bond dim." => (4, 64), "k" => 2) === true
    out = String(take!(io))
    t1 = now()
    @test count(==('\n'), out) == 1 && endswith(out, '\n')
    ts, body = split_line(chop(out))
    @test t0 <= ts <= t1                        # the local time of printing
    @test body == "CAMPS gate 3/10 (30%) bond dim.=4 (max 64) k=2"

    # No values: nothing after the percentage. Single values, floats and strings as `string`.
    log_progress!(p, 4)
    log_progress!(p, 5, "max S_E" => 0.5, "mode" => "xy", "n. of Paulis" => (120, 10^6))
    @test bodies(io) == ["CAMPS gate 4/10 (40%)",
                         "CAMPS gate 5/10 (50%) max S_E=0.5 mode=xy n. of Paulis=120 (max 1000000)"]

    # `indent` spaces between the timestamp and the text, as for the per-gate loops
    p = ProgressLog("PP gate", 4; io, frac = 0, indent = 2)
    log_progress!(p, 1, "n. of Paulis" => (7, 100))
    @test bodies(io) == ["  PP gate 1/4 (25%) n. of Paulis=7 (max 100)"]

    # log_message: same timestamp prefix, own indent
    log_message("p=0.020 CAMPS-PP sample 7/25"; io)
    log_message("stopped early"; io, indent = 3)
    @test bodies(io) == ["p=0.020 CAMPS-PP sample 7/25", "   stopped early"]
    @test log_message("x"; io = IOBuffer()) === nothing

    # Plain lines only: no carriage returns, ANSI escapes or cursor moves anywhere
    p = ProgressLog("x", 50; io, frac = 0)
    foreach(i -> log_progress!(p, i, "v" => (i, 50)), 1:50)
    log_done!(p, 50)
    out = String(take!(io))
    @test count(==('\n'), out) == 50
    @test !occursin('\r', out) && !occursin('\e', out)
  end

  @testset "percentage is floored" begin
    io = IOBuffer()
    p = ProgressLog("x", 1000; io, frac = 0)
    foreach(i -> log_progress!(p, i), (1, 9, 10, 995, 999, 1000))
    @test bodies(io) == ["x 1/1000 (0%)", "x 9/1000 (0%)", "x 10/1000 (1%)",
                         "x 995/1000 (99%)", "x 999/1000 (99%)", "x 1000/1000 (100%)"]
    p = ProgressLog("x", 3; io, frac = 0)
    foreach(i -> log_progress!(p, i), 1:3)
    @test bodies(io) == ["x 1/3 (33%)", "x 2/3 (66%)", "x 3/3 (100%)"]
  end

  @testset "10 % milestones" begin
    # Hand-derived: the first step at or past each multiple of n/10, then the last step.
    @test logged_steps(100, 1//10, 1:100) == 10:10:100
    @test logged_steps(30, 1//10, 1:30) == 3:3:30      # 30·0.1 = 3.0000000000000004 in Float64
    @test logged_steps(25, 1//10, 1:25) == [3, 5, 8, 10, 13, 15, 18, 20, 23, 25]
    @test logged_steps(7, 1//10, 1:7) == 1:7           # n < 10: every step passes a multiple
    @test logged_steps(3105, 1//10, 1:3105) ==         # gates per sample at N = 46, t = 23
          [311, 621, 932, 1242, 1553, 1863, 2174, 2484, 2795, 3105]
    @test logged_steps(1, 1//10, 1:1) == [1]
    @test logged_steps(10, 1//10, 0:10) == 1:10        # step 0 is not a milestone
    @test logged_steps(8, 1//4, 1:8) == [2, 4, 6, 8]
    @test logged_steps(10, 1//1, 1:10) == [10]         # frac ≥ 1: only the last step
    @test logged_steps(10, 2//1, 1:10) == [10]
    # Calls that skip steps: one line however many multiples were passed at once
    @test logged_steps(100, 1//10, [5, 35, 38, 40, 99, 100]) == [35, 40, 99, 100]
    @test logged_steps(100, 1//10, [10, 10, 10, 20]) == [10, 20]   # a step prints once
    # frac = 0 prints every (new) step
    @test logged_steps(5, 0//1, 1:5) == 1:5
    @test logged_steps(5, 0//1, [2, 2, 3]) == [2, 3]

    # Against the integer oracle, every n up to 400, several fractions (Float64 boundaries)
    bad = []
    for frac in (1//10, 1//20, 1//4, 1//3, 3//10, 1//1), n in 1:400
      got = logged_steps(n, frac, 1:n)
      got == oracle_steps(n, frac, 1:n) || push!(bad, (n, frac))
    end
    @test isempty(bad)
    # … and on random increasing call sequences (loops that report only some steps)
    rng = Xoshiro(7)
    for _ in 1:300
      n = rand(rng, 1:5000)
      calls = sort(unique(rand(rng, 0:n, rand(rng, 1:40))))
      frac = rand(rng, (1//10, 1//4, 1//3, 0//1))
      logged_steps(n, frac, calls) == oracle_steps(n, frac, calls) || push!(bad, (n, frac, calls))
    end
    @test isempty(bad)

    # A tiny frac neither overflows nor stops printing (every step is a new milestone)
    @test logged_steps(5, 1e-30, 1:5) == 1:5
    io = IOBuffer()
    p = ProgressLog("x", 10^9; io, frac = 1e-300, dt = Inf)
    @test log_progress!(p, 1) && log_progress!(p, 2) && log_progress!(p, 10^9)
  end

  @testset "dt trigger" begin
    io = IOBuffer()
    # dt = 0: a line is due at every new step whatever the milestones
    p = ProgressLog("x", 1000; io, frac = 0.5, dt = 0)
    @test all(i -> log_progress!(p, i), 1:20)
    @test !log_progress!(p, 20)                   # still once per step
    take!(io)
    # dt > 0: due once dt seconds have passed since the last line (or the creation).
    # (log_progress! is compiled by now, so the first call below is fast.)
    p = ProgressLog("x", 1000; io, frac = 0.5, dt = 0.5)
    @test !isdue(p, 1) && !log_progress!(p, 1)
    sleep(0.6)
    @test isdue(p, 2)
    @test log_progress!(p, 2, "k" => 3)
    @test !log_progress!(p, 3)                    # the line just printed restarts the clock
    @test log_progress!(p, 500)                   # milestones still print in between
    @test bodies(io) == ["x 2/1000 (0%) k=3", "x 500/1000 (50%)"]
  end

  @testset "log_done!" begin
    io = IOBuffer()
    # Loop stopped early (e.g. CAMPS reaching χ): the final line shows where it stopped
    p = ProgressLog("CAMPS gate", 100; io, dt = Inf)
    foreach(i -> log_progress!(p, i), 1:37)
    @test needs_final(p, 37)
    @test log_done!(p, 37, "bond dim." => (64, 64)) === nothing
    @test bodies(io) == ["CAMPS gate 10/100 (10%)", "CAMPS gate 20/100 (20%)",
                         "CAMPS gate 30/100 (30%)", "CAMPS gate 37/100 (37%) bond dim.=64 (max 64)"]
    @test !needs_final(p, 37)
    log_done!(p, 37)                              # a second call prints nothing
    @test isempty(bodies(io))

    # Loop that ran to the end: the last step was printed by log_progress!, no repeat
    p = ProgressLog("x", 20; io, dt = Inf)
    foreach(i -> log_progress!(p, i), 1:20)
    @test !needs_final(p, 20)
    log_done!(p, 20, "k" => 1)
    @test printed_step.(bodies(io)) == 2:2:20

    # Stopped exactly at a printed milestone: nothing more
    p = ProgressLog("x", 100; io, dt = Inf)
    foreach(i -> log_progress!(p, i), 1:30)
    log_done!(p, 30)
    @test printed_step.(bodies(io)) == [10, 20, 30]

    # A loop that applied nothing (stopped at step 0) still gets its line
    p = ProgressLog("x", 100; io, dt = Inf)
    log_done!(p, 0, "bond dim." => (2, 2))
    @test bodies(io) == ["x 0/100 (0%) bond dim.=2 (max 2)"]
  end

  @testset "disabled log" begin
    io = IOBuffer()
    p = ProgressLog("x", 10; io, frac = 0, dt = 0, enabled = false)
    @test !any(i -> isdue(p, i), 0:10)
    @test !any(i -> log_progress!(p, i, "k" => i), 1:10)
    @test !needs_final(p, 10) && !needs_final(p, 3)
    @test log_done!(p, 3, "k" => 1) === nothing
    @test isempty(take!(io))
  end

  @testset "n = 0" begin
    io = IOBuffer()
    p = ProgressLog("CAMPS gate", 0; io)          # e.g. `evolve` with t = 0
    @test isdue(p, 0) && needs_final(p, 0)
    @test log_progress!(p, 0, "k" => 0)
    log_done!(p, 0, "k" => 0)
    @test bodies(io) == ["CAMPS gate 0/0 (100%) k=0"]
    p = ProgressLog("CAMPS gate", 0; io)
    log_done!(p, 0)
    @test bodies(io) == ["CAMPS gate 0/0 (100%)"]
  end

  @testset "isdue matches log_progress! and has no side effects" begin
    for (n, frac) in ((100, 0.1), (25, 0.1), (7, 0.25), (40, 0.0))
      io = IOBuffer()
      p = ProgressLog("x", n; io, frac, dt = Inf)
      ok = true
      for i in 1:n
        d = isdue(p, i)
        before = position(io)
        ok &= isdue(p, i) == d                    # same answer twice, nothing printed
        ok &= position(io) == before
        ok &= log_progress!(p, i) == d
        ok &= !isdue(p, i)                        # a printed (or skipped) step is not due again
      end
      @test ok
    end
  end

  @testset "lines are flushed" begin
    mktempdir() do dir
      path = joinpath(dir, "log.txt")
      open(path, "w") do f
        # Control: a file IOStream buffers, so an unflushed line is not in the file yet
        println(f, "unflushed")
        @test filesize(path) == 0
        p = ProgressLog("x", 2; io = f, frac = 0)
        log_progress!(p, 1)
        @test countlines(path) == 2               # visible before close (the log flushed)
        log_message("m"; io = f)
        @test countlines(path) == 3
        log_done!(p, 2)
        @test countlines(path) == 4
      end
    end
  end

  @testset "constructor" begin
    p = ProgressLog("x", 10)
    @test p.io === stdout && p.frac == 0.1 && p.dt == 30.0 && p.enabled && p.indent == 0
    @test_throws ArgumentError ProgressLog("x", -1)
    @test_throws ArgumentError ProgressLog("x", 10; frac = -0.1)
    @test_throws ArgumentError ProgressLog("x", 10; dt = -1)
    @test_throws ArgumentError ProgressLog("x", 10; indent = -1)
  end

  # -- The evolution loops with showprogress = true --------------------------------------
  # Same results as with showprogress = false (bit for bit), nothing printed when it is
  # false, and the printed steps follow the milestone rule over the steps each loop runs.

  @testset "evolution loops print and give the same results" begin
    N, t, ones = 6, 4, [2, 5]
    gates, phases, _ = doped_xxz(N, t)
    M = length(gates)
    le = spec_layer_ends(N, t)
    O = test_observable(N, ones)
    ψ0 = computationalcamps(mt.QubitCAMPS, N, ones)

    function both(f)
      quiet, qlines = capture_stdout(() -> f(false))
      loud, llines = capture_stdout(() -> f(true))
      return quiet, loud, qlines, llines
    end

    @testset "camps_circuit_dynamics" begin
      χ = 4                                       # stops before the end
      q, l, ql, ll = both(sp -> camps_circuit_dynamics(ψ0, gates, phases, χ, O;
                                                       layer_ends = le, showprogress = sp, track = true))
      @test isempty(ql)
      @test isequal(q[2:end], l[2:end])           # k, stop gate, values, bond dimensions
      s = q[3]
      @test 0 < s < M
      @test check_loop_lines(ll, "CAMPS gate", M, ["bond dim.", "k"], 1:s)
      @test endswith(last(ll), " $s/$M ($(fld(100s, M))%) bond dim.=$(bonddim(q[1])) (max $χ) k=$(q[2]) (max $N)")
    end

    @testset "campspp_circuit_dynamics (CAMPS then PP)" begin
      χ = 4
      q, l, ql, ll = both(sp -> campspp_circuit_dynamics(ψ0, χ, 0.0, typemax(Int), gates, phases, O;
                                                         layer_ends = le, showprogress = sp, track = true))
      @test isempty(ql)
      @test isequal(q, l)
      s_switch = q[7]
      @test 0 < s_switch < M                      # both loops run
      camps = filter(x -> occursin("CAMPS gate", x), ll)
      ppl = filter(x -> occursin("PP gate", x), ll)
      @test length(camps) + length(ppl) == length(ll)
      @test check_loop_lines(camps, "CAMPS gate", M, ["bond dim.", "k"], 1:s_switch)
      @test check_loop_lines(ppl, "PP gate", M, ["n. of Paulis"], s_switch+1:M)
      @test findlast(x -> occursin("CAMPS gate", x), ll) < findfirst(x -> occursin("PP gate", x), ll)
    end

    @testset "pauliprop_circuit_dynamics (PP alone)" begin
      q, l, ql, ll = both(sp -> pauliprop_circuit_dynamics(ones, gates, phases, 0.0, typemax(Int), O;
                                                           layer_ends = le, showprogress = sp, track = true))
      @test isempty(ql)
      @test isequal(q, l)
      @test check_loop_lines(ll, "PP gate", M, ["n. of Paulis"], 1:M)
    end

    @testset "campssrc_circuit_dynamics" begin
      q, l, ql, ll = both(sp -> campssrc_circuit_dynamics(ψ0, gates, phases, 1e-11, O;
                                                          layer_ends = le, showprogress = sp, track = true))
      @test isempty(ql)
      @test isequal(q[2:end], l[2:end])
      @test check_loop_lines(ll, "CAMPS (SVD) gate", M, ["bond dim."], 1:M)
    end

    @testset "mps_circuit_dynamics" begin
      # A cutoff that truncates, so a gauge change from the progress display (entropies)
      # would change the values
      thl_mps = 1e-3
      q, l, ql, ll = both(sp -> mps_circuit_dynamics(basis_mps(N, ones), gates, phases, thl_mps, O, nothing;
                                                      layer_ends = le, showprogress = sp))
      @test isempty(ql)
      @test isequal(q[2], l[2])                   # recorded values
      @test check_loop_lines(ll, "MPS gate", M, ["bond dim.", "max S_E"], 1:M)
    end

    @testset "evolve, evolve_bonddim, evolve_deepcliffords" begin
      paulis = [PauliOperator(getpauli(g, N)) for g in gates]
      ψ = computationalcamps(mt.QubitCAMPS, N, ones)
      q, l, ql, ll = both(sp -> evolve(deepcopy(ψ), M, paulis, phases; showprogress = sp))
      @test isempty(ql) && q[2] == l[2]
      @test isequal(state_expectation(q[1], O), state_expectation(l[1], O))
      @test check_loop_lines(ll, "CAMPS gate", M, ["k"], 1:M)

      q, l, ql, ll = both(sp -> evolve(deepcopy(ψ), 0, paulis, phases; showprogress = sp))
      @test isempty(ql) && check_loop_lines(ll, "CAMPS gate", 0, ["k"], 0:0)

      q, l, ql, ll = both(sp -> DisentangleCAMPS.evolve_bonddim(deepcopy(ψ), 4, paulis, phases; showprogress = sp))
      @test isempty(ql) && q[2:3] == l[2:3]
      @test 0 < q[3] < M
      @test check_loop_lines(ll, "CAMPS gate", M, ["bond dim."], 1:q[3])

      tl = 5
      q, l, ql, ll = both(sp -> evolve_deepcliffords(Xoshiro(3), deepcopy(ψ), tl, paulis, phases; showprogress = sp))
      @test isempty(ql) && q[2] == l[2]
      @test isequal(state_expectation(q[1], O), state_expectation(l[1], O))
      @test check_loop_lines(ll, "CAMPS layer", tl, ["k"], 1:tl)
    end
  end
end
