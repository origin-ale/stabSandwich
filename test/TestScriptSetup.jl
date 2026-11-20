# Tests for scripts/ScriptSetup.jl, the shared setup of the top-level scripts (Phase 7). It
# lives outside CampsPP, so it is tested here in the root environment:
#   julia --project=. -t 2 test/TestScriptSetup.jl   (also run by test/runtests.jl)
#
# Spec (plan, Phase 7 decisions, and the module's docstrings): `key=value` overrides for
# every parameter of a script's `@parameters` block, parsed to the default's type (lists
# comma-separated, ranges a:b / a:s:b, Symbols with or without a leading `:`), applied right
# after each assignment so later lines see them; unknown, repeated or malformed keys are
# errors; `help` lists the parameters (as typed on the command line) and exits.
# `param_header` gives the `"name" => "value"` string pairs written into the file headers.
#
# In-process tests set `ARGS` around each block. `help` calls `exit(0)` outside an
# interactive session, so it (and the plain script path) is tested in child processes.
# Untested: the interactive branch of `help` (throws instead of exiting).

using Test
using LinearAlgebra: BLAS
import Strided

const ROOT = dirname(@__DIR__)
const SETUP = joinpath(ROOT, "scripts", "ScriptSetup.jl")
include(SETUP)
using .ScriptSetup

"The exception inside the `LoadError`s that `include_string` wraps it in."
unwrap(e) = e isa LoadError ? unwrap(e.error) : e

"""
The globals of a module, read in the latest world: they are created by `include_string`
after the calling function started, so a plain `m.x` would not see them (Julia ≥ 1.12).
"""
struct Globals
  m::Module
end
Base.getproperty(g::Globals, s::Symbol) = Base.invokelatest(getglobal, getfield(g, :m), s)

"""
Evaluate `@parameters begin <block> end` in a fresh module with `ARGS = args`, as a script
would at top level. Returns the module's globals (the parameters, as `Globals`) and the
returned pairs. Errors are rethrown unwrapped.
"""
function run_params(block::AbstractString, args::Vector{String} = String[])
  saved = copy(ARGS)
  m = Module(:ParamsTest)
  Core.eval(m, :(const ScriptSetup = $ScriptSetup))
  Core.eval(m, :(using .ScriptSetup))
  try
    empty!(ARGS); append!(ARGS, args)
    params = include_string(m, "@parameters begin\n$block\nend")
    return Globals(m), params
  catch e
    throw(unwrap(e))
  finally
    empty!(ARGS); append!(ARGS, saved)
  end
end

"Value of the single parameter `x = <default>` after the override `x=<arg>`."
parsed(default::AbstractString, arg::AbstractString) = run_params("x = $default", ["x=$arg"])[1].x

"Whether `f()` throws an `ErrorException` whose message contains every one of `parts`."
function throws_error(f, parts...)
  try
    f()
  catch e
    e isa ErrorException || return false
    return all(p -> occursin(p, e.msg), parts)
  end
  return false
end

"Whether `f()` throws an `ArgumentError` (unwrapped from `LoadError`s) whose message contains `part`."
function unwrap_arg(f, part)
  try
    f()
  catch e
    e = unwrap(e)
    return e isa ArgumentError && occursin(part, e.msg)
  end
  return false
end

@testset "ScriptSetup.jl" begin
  @testset "parsing per type" begin
    # Int: underscores, integral floats and expressions; non-integers are errors
    @test parsed("46", "12") === 12
    @test parsed("46", "-3") === -3
    @test parsed("46", "1_000") === 1000
    @test parsed("46", " 12 ") === 12
    @test parsed("10", "1e7") === 10_000_000
    @test parsed("10", "2^10") === 1024
    @test parsed("10", "3*4") === 12
    for bad in ("1.5", "abc", "", "1e-3", "sin(1)", "1+")
      @test throws_error(() -> parsed("46", bad), "Cannot parse `x=$bad`", "Int64")
    end

    # Float64: literals, π/pi and + - * / ^ of them (as Julia evaluates the same text)
    @test parsed("0.6", "0.04") === 0.04
    @test parsed("0.6", "2") === 2.0
    @test parsed("0.6", "-0.5") === -0.5
    @test parsed("1e-10", "1e-12") === 1e-12
    @test parsed("π/4", "pi/4") === π/4
    @test parsed("π/4", "3pi/16") === 3π/16
    @test parsed("π/4", "3π/16") === 3π/16
    @test parsed("π/4", "-pi/3") === -π/3
    @test parsed("π/4", "pi") === Float64(π)
    for bad in ("x", "sin(1)", "pi/", "", "0.1.2", "1,2")
      @test throws_error(() -> parsed("0.6", bad), "Cannot parse `x=$bad`", "Float64")
    end
    # An Irrational default (`ϕ = π`) takes a Float64 override; without one it stays π
    @test parsed("π", "pi/2") === π/2
    @test run_params("ϕ = π")[1].ϕ === π

    # Bool
    for (s, v) in ("true" => true, "TRUE" => true, "1" => true, "yes" => true,
                   "false" => false, "0" => false, "no" => false, "False" => false)
      @test parsed("false", s) === v
    end
    @test throws_error(() -> parsed("false", "maybe"), "Cannot parse `x=maybe`", "Bool")

    # Symbol: with or without the colon
    @test parsed(":chi3", "entangle") === :entangle
    @test parsed(":chi3", ":entangle") === :entangle
    @test parsed(":xy", "z") === :z

    # String: kept verbatim, including `=` after the first one and an empty value
    @test parsed("\"output/X\"", "/tmp/run 1") == "/tmp/run 1"
    @test parsed("\"a\"", "k=v=w") == "k=v=w"
    @test parsed("\"a\"", "") == ""
    @test parsed("\"a\"", "12") == "12"

    # Vectors: comma-separated, element type of the default, empty value → empty vector
    v = parsed("[0.0, 0.01]", "0,0.04")
    @test v == [0.0, 0.04] && v isa Vector{Float64}
    @test parsed("[0.0]", "pi/4,3pi/16") == [π/4, 3π/16]
    v = parsed("[:campspp]", "campspp,pp")
    @test v == [:campspp, :pp] && v isa Vector{Symbol}
    @test parsed("[:campspp]", ":camps,:mps") == [:camps, :mps]
    v = parsed("[12, 24]", "8")
    @test v == [8] && v isa Vector{Int}
    v = parsed("[0.0, 0.01]", "")
    @test isempty(v) && v isa Vector{Float64}
    @test parsed("[\"a\"]", "x,y") == ["x", "y"]
    @test throws_error(() -> parsed("[1, 2]", "1,,2"), "Cannot parse `x=1,,2`")
    @test throws_error(() -> parsed("[1, 2]", "1.5,2"), "Cannot parse `x=1.5,2`")
    @test throws_error(() -> parsed("[]", "1"), "Cannot parse `x=1`", "Vector{Any}")

    # Ranges: start:stop or start:step:stop
    @test parsed("4:2:8", "4:4:16") === 4:4:16
    @test parsed("4:2:8", "12:24") == 12:24
    @test parsed("1:3", "2:5") === 2:5
    @test throws_error(() -> parsed("4:2:8", "4"), "Cannot parse `x=4`", "range")
    @test throws_error(() -> parsed("4:2:8", "1:2:3:4"), "Cannot parse `x=1:2:3:4`")

    # Types that cannot come from the command line: an error only when overridden
    @test throws_error(() -> parsed("(1, 2)", "3"), "Cannot parse `x=3`", "cannot be set from the command line")
    @test run_params("x = (1, 2)")[1].x == (1, 2)
  end

  @testset "fixed parser bugs (Phase 7 test findings, fixed 2026-11-18)" begin
    # 1. Negative powers for a float parameter: Float64 arithmetic, not Int powers (was a
    #    raw DomainError). Large integer literals and powers become floats, not wrapped Ints.
    @test parsed("0.6", "2^-1") === 0.5
    @test parsed("0.6", "1e-12") === 1e-12
    @test parsed("0.6", "-2^-2") === -0.25
    @test parsed("0.6", "10^20") === 1e20
    @test parsed("0.6", "9223372036854775808") === 9.223372036854776e18        # Int128 literal
    @test parsed("0.6", "1" * "0"^41) === 1e41                                 # BigInt literal
    @test parsed("[0.0]", "2^-1,pi^2") == [0.5, π^2]
    #    Any remaining math error is reported as a parse error naming the key, not raw
    @test throws_error(() -> parsed("0.6", "(-8)^(1/3)"), "Cannot parse `x=(-8)^(1/3)`", "DomainError")

    # 2. Integer parameters: exact up to the Int64 limits, an error beyond (was a silent
    #    wrap-around: 2^100 → 0, 10^19 → negative)
    @test parsed("46", "2^62") === 2^62
    @test parsed("46", "2^63-1") === typemax(Int)
    @test parsed("46", "-2^63") === typemin(Int)
    @test parsed("46", "9223372036854775807") === typemax(Int)
    @test parsed("46", "4/2") === 2
    for bad in ("2^63", "2^100", "10^19", "-2^63-1", "9223372036854775808", "-9223372036854775809",
                "1" * "0"^41, "3*10^18*4")
      @test throws_error(() -> parsed("46", bad), "Cannot parse `x=$bad` as Int64", "out of range for Int64")
    end
    @test throws_error(() -> parsed("10^7", "1e19"), "Cannot parse `x=1e19`", "out of range")
    @test throws_error(() -> parsed("46", "2^-1"), "Cannot parse `x=2^-1`", "not an integer")
    @test throws_error(() -> parsed("46", "0^-1"), "Cannot parse `x=0^-1`", "not an integer")
    @test throws_error(() -> parsed("46", "1/0"), "Cannot parse `x=1/0`", "not an integer")
    @test throws_error(() -> parsed("[1]", "1,2^64"), "Cannot parse `x=1,2^64`", "out of range")

    # 3. Spaces around list items are dropped for every element type (a quoted
    #    `methods="campspp, pp"`; was Symbol(" pp")), and around a Symbol value
    v = parsed("[:campspp]", "campspp, pp")
    @test v == [:campspp, :pp] && v isa Vector{Symbol}
    @test parsed("[:campspp]", " :camps , :mps ") == [:camps, :mps]
    @test parsed("[\"a\"]", "x, y") == ["x", "y"]
    @test parsed("[0.0]", "0, pi/4") == [0.0, π/4]
    @test parsed("[1]", " 12 , 24") == [12, 24]
    @test parsed(":chi3", " entangle ") === :entangle
    @test parsed(":chi3", " :entangle") === :entangle
    #    A scalar String is still kept verbatim, and an empty number item is still an error
    @test parsed("\"a\"", " a b ") == " a b "
    @test throws_error(() -> parsed("[1, 2]", "1, ,2"), "Cannot parse `x=1, ,2`", "empty number")
  end

  @testset "exact number expressions and empty items (decided 2026-11-18)" begin
    # Integer literals are exact rationals through + - * / ^, converted to the parameter's
    # type once at the end, so a float result is the correctly rounded value of the exact
    # expression (was Float64 `^`: 10^-12 came out 1 ulp above 1e-12)
    for k in 1:20
      @test parsed("0.6", "10^-$k") === parse(Float64, "1e-$k")
    end
    @test parsed("0.6", "1/3") === 1/3
    @test parsed("0.6", "1/3+1/6") === 0.5                  # exact: 1/3 + 1/6 = 1/2
    @test parsed("0.6", "2^-1") === 0.5
    @test parsed("0.6", "(1/3)*3") === 1.0
    @test parsed("0.6", "2^-5000") === 0.0
    #    A float literal or π turns the operands into Float64 first, as in Julia source
    @test parsed("0.6", "3*0.1") === 3*0.1                  # 0.30000000000000004, not 0.3
    @test parsed("0.6", "0.1+0.2") === 0.1+0.2
    @test parsed("π/4", "3pi/16") === 3π/16
    @test parsed("π/4", "pi/4") === π/4
    @test parsed("0.6", "2^0.5") === sqrt(2.0)
    #    Division by zero as in Julia (±Inf), except 0/0, which is an error (Julia: NaN)
    @test parsed("0.6", "1/0") === Inf
    @test parsed("0.6", "-1/0") === -Inf
    @test parsed("0.6", "0^-1") === Inf
    @test throws_error(() -> parsed("0.6", "0/0"), "Cannot parse `x=0/0`")
    @test throws_error(() -> parsed("46", "0/0"), "Cannot parse `x=0/0`")
    #    Huge exponents (> 4096 in absolute value) give Float64 overflow, not a huge BigInt
    @test parsed("0.6", "2^5000") === Inf
    @test throws_error(() -> parsed("46", "2^5000"), "Cannot parse `x=2^5000`")

    # Integer parameters: an exact integral result is accepted, whatever the route
    @test parsed("46", "100/2") === 50
    @test parsed("46", "6/4*2") === 3
    @test parsed("46", "(2^63-1)/1") === typemax(Int)
    @test throws_error(() -> parsed("46", "2^62*2"), "Cannot parse `x=2^62*2`", "out of range for Int64")
    @test throws_error(() -> parsed("46", "1e30"), "Cannot parse `x=1e30`", "out of range for Int64")
    @test throws_error(() -> parsed("46", "3/2"), "Cannot parse `x=3/2`", "not an integer")

    # Empty Symbols and empty Symbol/String list items are errors (were Symbol("") / "")
    for bad in ("", ":", " ")
      @test throws_error(() -> parsed(":chi3", bad), "Cannot parse `x=$bad`", "empty Symbol")
    end
    for (default, bad) in (("[:a]", "a,,b"), ("[:a]", "a,"), ("[:a]", ",a"), ("[:a]", "a, ,b"),
                           ("[\"a\"]", "a,,b"), ("[\"a\"]", "a, "))
      @test throws_error(() -> parsed(default, bad), "Cannot parse `x=$bad`", "empty item in the list `$bad`")
    end
    @test throws_error(() -> parsed("[1]", "1,,2"), "Cannot parse `x=1,,2`", "empty number")
    #    Unchanged: an empty String scalar is "", an empty value is an empty list
    @test parsed("\"a\"", "") === ""
    v = parsed("[:a]", "")
    @test isempty(v) && v isa Vector{Symbol}
    v = parsed("[\"a\"]", "")
    @test isempty(v) && v isa Vector{String}
  end

  @testset "@parameters: order, globals, return value" begin
    block = """
      N = 46
      t = N ÷ 2                     # derived: sees an override of N
      μ = 0.6
      magic_probs = [0.0, 0.01]
      methods = [:campspp]
      sites = [1, N]
      outdir = "output/Script"
      """
    m, params = run_params(block)
    @test (m.N, m.t, m.μ, m.sites) == (46, 23, 0.6, [1, 46])
    @test params isa Vector{Pair{Symbol, Any}}
    @test first.(params) == [:N, :t, :μ, :magic_probs, :methods, :sites, :outdir]
    @test last.(params) == [46, 23, 0.6, [0.0, 0.01], [:campspp], [1, 46], "output/Script"]

    # An override applies right after its own line: later lines see it, earlier ones don't
    m, params = run_params(block, ["N=12"])
    @test (m.N, m.t, m.sites) == (12, 6, [1, 12])
    # A derived parameter can itself be overridden (the override wins over the formula)
    m, params = run_params(block, ["t=5", "N=12"])     # order on the command line irrelevant
    @test (m.N, m.t) == (12, 5)
    @test last.(params)[1:2] == [12, 5]
    # Overriding a later line does not change earlier ones
    m, _ = run_params("a = 1\nb = a + 1\nc = b * 10", ["b=7"])
    @test (m.a, m.b, m.c) == (1, 7, 70)

    # Several overrides of every kind at once; the returned pairs hold the overridden values
    m, params = run_params(block, ["methods=campspp,pp", "magic_probs=0,0.04", "μ=0.5", "outdir=/tmp/o"])
    @test Dict(params)[:methods] == [:campspp, :pp]
    @test Dict(params)[:magic_probs] == [0.0, 0.04]
    @test (m.μ, m.outdir) == (0.5, "/tmp/o")

    # The micro sign µ (U+00B5), as typed on many keyboards, names the parameter μ (U+03BC)
    m, _ = run_params("μ = 0.6", ["µ=0.3"])
    @test m.μ == 0.3

    # Comments, blank lines and line numbers in the block are fine; values are any expression
    m, params = run_params("# comment\n\nx = 2 + 3   # five\n\ny = string(:a)\n")
    @test params == [:x => 5, :y => "a"]
  end

  @testset "@parameters: errors" begin
    block = "N = 46\nt = N ÷ 2"
    @test throws_error(() -> run_params(block, ["Nx=3"]), "Unrecognised parameter(s): Nx.", "Available: N, t")
    @test throws_error(() -> run_params(block, ["b=1", "a=2", "N=4"]), "Unrecognised parameter(s): a, b.")
    @test throws_error(() -> run_params(block, ["N=3", "N=4"]), "Parameter `N` given twice")
    @test throws_error(() -> run_params("μ = 0.6", ["μ=0.5", "µ=0.4"]), "given twice")
    @test throws_error(() -> run_params(block, ["N"]), "Argument `N` is not of the form key=value")
    @test throws_error(() -> run_params(block, ["=3"]), "Argument `=3` is not of the form key=value")
    @test throws_error(() -> run_params(block, ["N=abc"]), "Cannot parse `N=abc` as Int64 (default 46)")

    # Malformed blocks are rejected when the macro expands
    @test unwrap_arg(() -> run_params("N = 1\nN = 2"), "`N` assigned twice")
    @test unwrap_arg(() -> run_params("f(x) = 1"), "only `name = value` lines")
    @test unwrap_arg(() -> run_params("N += 1"), "only `name = value` lines")
    @test unwrap_arg(() -> run_params("println(1)"), "only `name = value` lines")
    m = Module(:ParamsTest2)
    Core.eval(m, :(const ScriptSetup = $ScriptSetup)); Core.eval(m, :(using .ScriptSetup))
    @test unwrap_arg(() -> include_string(m, "@parameters N = 3"), "expects a begin … end block")
  end

  @testset "param_header" begin
    params = Pair{Symbol, Any}[:N => 46, :magic_probs => [0.0, 0.01], :methods => [:campspp, :pp],
      :camps_crit => :chi3, :outdir => "output/X", :Ns => 4:2:8, :ts => 1:3, :ϕ => π/4,
      :flag => true, :thl => 1e-10, :empty => Float64[]]
    h = param_header(params)
    @test h isa Vector{Pair{String, String}}
    @test h == ["N" => "46", "magic_probs" => "0.0,0.01", "methods" => "campspp,pp",
                "camps_crit" => "chi3", "outdir" => "output/X", "Ns" => "4:2:8", "ts" => "1:3",
                "ϕ" => "0.7853981633974483", "flag" => "true", "thl" => "1.0e-10", "empty" => ""]
    @test !any(p -> occursin(' ', p.second), h)        # each header entry stays one k=v word
    @test first.(param_header(params; exclude = (:outdir, :Ns))) ==
          ["N", "magic_probs", "methods", "camps_crit", "ts", "ϕ", "flag", "thl", "empty"]
    @test param_header(params; exclude = [:N])[1] == ("magic_probs" => "0.0,0.01")
    @test isempty(param_header(Pair{Symbol, Any}[]))
    # Straight from @parameters, with overrides
    _, params = run_params("N = 46\nmethods = [:campspp]\nouth = \"o\"", ["methods=pp,mps"])
    @test param_header(params; exclude = (:outh,)) == ["N" => "46", "methods" => "pp,mps"]
  end

  @testset "setup_threads" begin
    # The scripts call it with 1 or with their `blas_threads` parameter (default 1). Tested at
    # n = 1 and n = Threads.nthreads() (Strided caps its count at nthreads, so a larger n is
    # not tested). runtests.jl runs this file with `-t 2`.
    blas0, strided0 = BLAS.get_num_threads(), Strided.get_num_threads()
    try
      for n in (1, Threads.nthreads(), 1)       # back to 1: each call sets the count
        @test setup_threads(n) === n
        @test BLAS.get_num_threads() == n && Strided.get_num_threads() == n
      end
      n = BLAS.get_num_threads()
      @test_throws ArgumentError setup_threads(0)
      @test BLAS.get_num_threads() == n                  # unchanged by the rejected call
    finally
      BLAS.set_num_threads(blas0); Strided.set_num_threads(strided0)
    end
  end

  @testset "as a script (child processes)" begin
    mktempdir() do dir
      script = joinpath(dir, "MyScript.jl")
      write(script, """
        @isdefined(ScriptSetup) || include($(repr(SETUP)))
        using .ScriptSetup
        params = @parameters begin
          N = 46
          magic_probs = [0.0, 0.01]
          methods = [:campspp, :pp]
          camps_crit = :chi3
          outdir = "output/X"
          Ns = 4:2:8
        end
        println("RAN ", join(["\$k=\$v" for (k, v) in param_header(params)], " "))
        """)
      run_script(args...) = begin
        obuf, ebuf = IOBuffer(), IOBuffer()
        cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no $script $args`
        proc = run(pipeline(ignorestatus(cmd); stdout = obuf, stderr = ebuf))
        return proc.exitcode, String(take!(obuf)), String(take!(ebuf))
      end

      # Plain run with overrides from the real command line
      code, out, _ = run_script("N=12", "methods=pp", "camps_crit=:entangle")
      @test code == 0
      @test out == "RAN N=12 magic_probs=0.0,0.01 methods=pp camps_crit=entangle outdir=output/X Ns=4:2:8\n"
      # A quoted list with spaces (one argument) and a large integer expression
      code, out, _ = run_script("methods=pp, campspp", "magic_probs=0, 2^-5", "N=2^40")
      @test code == 0
      @test out == "RAN N=1099511627776 magic_probs=0.0,0.03125 methods=pp,campspp camps_crit=chi3 outdir=output/X Ns=4:2:8\n"
      code, out, err = run_script("N=2^64")
      @test code != 0 && !occursin("RAN", out) && occursin("out of range for Int64", err)

      # help: lists the values (after overrides) as typed on the command line, exits 0,
      # and the script body does not run
      code, out, _ = run_script("help", "N=12")
      @test code == 0
      @test out == """
        Parameters (override with key=value; lists comma-separated, ranges a:b or a:s:b):
          N           = 12
          magic_probs = 0.0,0.01
          methods     = campspp,pp
          camps_crit  = chi3
          outdir      = output/X
          Ns          = 4:2:8
        """
      @test run_script("--help")[1:2] == run_script("-h")[1:2] == (0, replace(out, "= 12" => "= 46"))

      # An unknown key is an error even with help (checked first); nonzero exit
      code, out, err = run_script("help", "Nx=1")
      @test code != 0 && !occursin("RAN", out) && occursin("Unrecognised parameter(s): Nx", err)
    end
  end
end
