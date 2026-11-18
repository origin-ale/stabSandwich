"""
Shared setup for the top-level simulation scripts (not part of CampsPP, which stays a
standalone package). Use from a script in the repository root:

    @isdefined(ScriptSetup) || include(joinpath(@__DIR__, "scripts", "ScriptSetup.jl"))
    using .ScriptSetup

(the guard lets several scripts be `include`d in one session without redefining the module)

    params = @parameters begin
      N = 46                      # system size
      magic_probs = [0.0, 0.01]   # doping probabilities
      methods = [:campspp]
      outdir = "output/MyScript"
    end
    setup_threads(Threads.nthreads())
    initialize_output(joinpath(outdir, "x.txt"), "obs", param_header(params))

and run with `key=value` overrides, e.g.

    julia --project=. MyScript.jl N=12 magic_probs=0,0.04 methods=campspp,pp

`julia --project=. MyScript.jl help` lists the parameters and their values.
Symbols are given without the colon (`camps_crit=entangle`; a leading `:` is also accepted).
"""
module ScriptSetup

using LinearAlgebra: BLAS
import Strided

export @parameters, param_header, setup_threads

# == Threads ============================================================================

"""```setup_threads(n)```

Give BLAS and Strided (the ITensors tensor backend) ```n``` threads each. This is what \
the scripts did by hand before (```Strided.disable_threads()``` followed by \
```set_num_threads(n)```, whose net effect is ```n```). Returns ```n```."""
function setup_threads(n::Integer)
  n >= 1 || throw(ArgumentError("need at least 1 thread, got $n"))
  BLAS.set_num_threads(n)
  Strided.set_num_threads(n)
  return n
end

# == Command-line parameters ============================================================

const HELP_WORDS = ("help", "-h", "--help")

# Keys are normalized as Julia normalizes identifiers (NFC plus its character mappings, e.g.
# the micro sign µ U+00B5 → μ U+03BC), by parsing them as one; so they match the variable
# names in the source. A key that is not an identifier is kept as is (and then unknown).
function normalize_key(key::AbstractString)
  ex = Meta.parse(key; raise = false)
  return ex isa Symbol ? string(ex) : String(key)
end

"""Overrides given on the command line: `key => (raw string, used flag)`."""
mutable struct Overrides
  values::Dict{String, String}
  used::Set{String}
  help::Bool
end

function parse_args(args)
  values = Dict{String, String}()
  help = false
  for arg in args
    if arg in HELP_WORDS
      help = true
      continue
    end
    kv = split(arg, "="; limit = 2)
    (length(kv) == 2 && !isempty(kv[1])) ||
      error("Argument `$arg` is not of the form key=value")
    key = normalize_key(kv[1])
    haskey(values, key) && error("Parameter `$(kv[1])` given twice")
    values[key] = kv[2]
  end
  return Overrides(values, Set{String}(), help)
end

function override(ov::Overrides, name::Symbol, default)
  key = normalize_key(string(name))
  haskey(ov.values, key) || return default
  push!(ov.used, key)
  str = ov.values[key]
  try
    return parse_value(typeof(default), str)
  catch e
    e isa ArgumentError || rethrow()
    error("Cannot parse `$key=$str` as $(typeof(default)) (default $(repr(default))): " *
          sprint(showerror, e))
  end
end

function check_unused(ov::Overrides, names)
  unknown = setdiff(keys(ov.values), ov.used)
  isempty(unknown) && return nothing
  error("Unrecognised parameter(s): $(join(sort(collect(unknown)), ", ")). " *
        "Available: $(join(names, ", "))")
end

# Numbers: literals, `π`/`pi` and + - * / ^ of those (e.g. `pi/4`, `3pi/16`, `1e-10`).
const NUM_OPS = Dict(:+ => +, :- => -, :* => *, :/ => /, :^ => ^)

evalnum(x::Number) = x
function evalnum(x::Symbol)
  x in (:π, :pi) && return π
  throw(ArgumentError("unknown name `$x` in a number"))
end
function evalnum(ex::Expr)
  if ex.head == :call && haskey(NUM_OPS, ex.args[1])
    return NUM_OPS[ex.args[1]](map(evalnum, ex.args[2:end])...)
  end
  throw(ArgumentError("`$ex` is not a number expression"))
end
evalnum(x) = throw(ArgumentError("`$x` is not a number"))

function parsenum(::Type{T}, str::AbstractString) where {T <: Real}
  s = replace(strip(str), "_" => "")
  isempty(s) && throw(ArgumentError("empty number"))
  T <: Integer && (n = tryparse(T, s)) !== nothing && return n
  ex = try
    Meta.parse(s)
  catch
    throw(ArgumentError("`$s` is not a number"))
  end
  x = evalnum(ex)
  # Integers also from an integral float or expression, e.g. 1e7 or 2^10
  if T <: Integer && !(x isa Integer)
    isinteger(x) || throw(ArgumentError("`$s` is not an integer"))
  end
  return convert(T, x)
end

parse_value(::Type{T}, str::AbstractString) where {T <: AbstractString} = String(str)
parse_value(::Type{Symbol}, str::AbstractString) = Symbol(chopprefix(str, ":"))
function parse_value(::Type{Bool}, str::AbstractString)
  s = lowercase(strip(str))
  s in ("true", "1", "yes") && return true
  s in ("false", "0", "no") && return false
  throw(ArgumentError("`$str` is not a Bool"))
end
parse_value(::Type{T}, str::AbstractString) where {T <: Real} = parsenum(T, str)
# A real default such as `ϕ = π/4` is an Irrational only when written as `π` alone
parse_value(::Type{<:AbstractIrrational}, str::AbstractString) = parsenum(Float64, str)

# Lists: comma-separated, no spaces; an empty value gives an empty list
function parse_value(::Type{V}, str::AbstractString) where {V <: AbstractVector}
  T = eltype(V)
  T === Any && throw(ArgumentError("cannot infer the element type of a Vector{Any} default"))
  isempty(str) && return T[]
  return T[parse_value(T, s) for s in split(str, ",")]
end

# Ranges: start:stop or start:step:stop
function parse_value(::Type{R}, str::AbstractString) where {R <: AbstractRange}
  T = eltype(R)
  parts = [parse_value(T, s) for s in split(str, ":")]
  length(parts) == 2 && return parts[1]:parts[2]
  length(parts) == 3 && return parts[1]:parts[2]:parts[3]
  throw(ArgumentError("`$str` is not a range start:stop or start:step:stop"))
end

parse_value(::Type{T}, str::AbstractString) where {T} =
  throw(ArgumentError("parameters of type $T cannot be set from the command line"))

# Written as it would be typed on the command line
help_value(v) = repr(v)
help_value(v::Symbol) = string(v)
help_value(v::AbstractString) = string(v)
help_value(v::AbstractVector) = join(map(help_value, v), ",")
help_value(v::AbstractRange) = string(v)

function print_help(io::IO, params)
  println(io, "Parameters (override with key=value; lists comma-separated, ranges a:b or a:s:b):")
  w = maximum(p -> textwidth(string(p.first)), params; init = 0)
  for (k, v) in params
    println(io, "  ", rpad(string(k), w), " = ", help_value(v))
  end
end

function help_and_stop(params)
  print_help(stdout, params)
  isinteractive() && error("@parameters: `help` given, parameters listed above; script not run")
  exit(0)
end

"""```params = @parameters begin … end```

Evaluate a block of ```name = value``` assignments as global variables, in order, \
replacing each value by its ```name=value``` command-line override (from ```ARGS```) \
right after its assignment, so later lines see overridden values (e.g. ```t = N ÷ 2``` \
follows ```N=…```). Overrides are parsed to the default's type: numbers (also \
```pi/4```-style expressions), ```Bool```, ```Symbol```, ```String```, comma-separated \
vectors, ```a:b```/```a:s:b``` ranges. Unknown keys are an error. With ```help``` on \
the command line, prints the parameters and their values and exits (in an interactive \
session it throws instead, so the session stays open). Returns the \
parameters as a vector of ```name => value``` pairs in block order (see \
```param_header```). Only plain assignments are allowed in the block."""
macro parameters(block)
  (block isa Expr && block.head == :block) ||
    throw(ArgumentError("@parameters expects a begin … end block"))
  ov = gensym("overrides")
  names = Symbol[]
  body = Expr(:block, :($ov = $parse_args(Base.ARGS)))
  for ex in block.args
    if ex isa LineNumberNode
      push!(body.args, ex)
    elseif ex isa Expr && ex.head == :(=) && ex.args[1] isa Symbol
      name = ex.args[1]
      name in names && throw(ArgumentError("@parameters: `$name` assigned twice"))
      push!(names, name)
      push!(body.args, ex, :($name = $override($ov, $(QuoteNode(name)), $name)))
    else
      throw(ArgumentError("@parameters: only `name = value` lines are allowed, got `$ex`"))
    end
  end
  params = gensym("params")
  push!(body.args,
    :($check_unused($ov, $names)),
    :($params = Pair{Symbol, Any}[$([:($(QuoteNode(n)) => $n) for n in names]...)]),
    :($ov.help && $help_and_stop($params)),
    params)
  return esc(body)
end

header_value(v) = string(v)
header_value(v::AbstractVector) = join(map(header_value, v), ",")
header_value(v::AbstractRange) = string(v)

"""```param_header(params; exclude = ())```

The parameters returned by ```@parameters``` as ```"name" => "value"``` string pairs \
for ```initialize_output```, in block order, without the names in ```exclude```. \
Vectors are written comma-separated without spaces, so that each header entry stays a \
single ```k=v``` word."""
param_header(params; exclude = ()) =
  [string(k) => header_value(v) for (k, v) in params if !(k in exclude)]

end
