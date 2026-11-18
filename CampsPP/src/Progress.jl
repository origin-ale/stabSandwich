using Dates

# == Plain-text progress output =========================================================
#
# Complete, timestamped, flushed lines: no `\r`, ANSI codes or cursor moves, so a log
# redirected to a file (e.g. by scripts/remote/run.sh) can be followed with `tail -f`.

timestamp() = Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS")

"""```log_message(msg; io = stdout, indent = 0)```

Print ```[yyyy-mm-dd HH:MM:SS] msg``` as one complete line (```indent``` spaces after the \
timestamp) and flush ```io```."""
function log_message(msg::AbstractString; io::IO = stdout, indent::Integer = 0)
  println(io, "[", timestamp(), "] ", " "^indent, msg)
  flush(io)
  return nothing
end

"""```ProgressLog(desc, n; io = stdout, frac = 0.1, dt = 30.0, enabled = true, indent = 0)```

Plain-text progress for a loop of at most ```n``` steps. Each line reads \
```[yyyy-mm-dd HH:MM:SS] desc i/n (pct%) key=value …```.

```log_progress!(p, i, values...)``` prints a line when one is due (see ```isdue```): \
at the last step, when the fraction done reaches the next multiple of ```frac```, or \
when ```dt``` seconds have passed since the last line. ```frac = 0``` prints every step. \
```log_done!(p, i, values...)``` prints the final line unless step ```i``` was already \
printed (e.g. a loop that stopped before ```n```).

Each value is a ```key => v``` pair, printed as ```key=v```; a tuple ```key => (v, vmax)``` \
is printed as ```key=v (max vmax)```. A disabled log prints nothing. Lines go to ```io``` \
(default ```stdout```, the same stream as the scripts' other messages)."""
mutable struct ProgressLog
  io::IO
  desc::String
  n::Int
  frac::Float64
  dt::Float64
  enabled::Bool
  indent::Int
  printed::Int          # last step printed (-1: none yet)
  milestone::Float64    # multiples of frac reached at the last printed step
  last_time::Float64    # time() of the last line, or of the creation
end

function ProgressLog(desc::AbstractString, n::Integer; io::IO = stdout, frac::Real = 0.1,
                     dt::Real = 30.0, enabled::Bool = true, indent::Integer = 0)
  n >= 0 || throw(ArgumentError("n must be ≥ 0, got $n"))
  frac >= 0 || throw(ArgumentError("frac must be ≥ 0, got $frac"))
  dt >= 0 || throw(ArgumentError("dt must be ≥ 0, got $dt"))
  indent >= 0 || throw(ArgumentError("indent must be ≥ 0, got $indent"))
  return ProgressLog(io, String(desc), n, frac, dt, enabled, indent, -1, 0.0, time())
end

# Number of multiples of p.frac in the fraction done i/n; the tolerance keeps exact multiples
# (e.g. 30/100 at frac 0.1, where 0.3/0.1 = 2.9999…) from landing one step late. A Float64,
# so that a tiny frac cannot overflow an Int.
milestone(p::ProgressLog, i::Integer) = floor(i / (p.n * p.frac) + 1e-9)

"""```isdue(p::ProgressLog, i)```

Whether ```log_progress!(p, i, …)``` would print a line at step ```i```. Use it to skip \
computing an expensive value when no line is due."""
function isdue(p::ProgressLog, i::Integer)
  p.enabled || return false
  i == p.printed && return false
  (i >= p.n || p.frac == 0) && return true
  milestone(p, i) > p.milestone && return true
  return time() - p.last_time >= p.dt
end

format_value(v) = string(v)
format_value(v::Tuple{Any, Any}) = string(v[1], " (max ", v[2], ")")

function print_line(p::ProgressLog, i::Integer, values)
  # floor, so that 100% means done (995/1000 shows 99%)
  pct = p.n == 0 ? 100 : floor(Int, 100 * i / p.n)
  buf = IOBuffer()
  print(buf, p.desc, " ", i, "/", p.n, " (", pct, "%)")
  for (k, v) in values
    print(buf, " ", k, "=", format_value(v))
  end
  log_message(String(take!(buf)); io = p.io, indent = p.indent)
  p.printed = i
  p.last_time = time()
  if p.frac > 0 && p.n > 0
    p.milestone = max(p.milestone, milestone(p, i))
  end
  return nothing
end

"""```log_progress!(p::ProgressLog, i, values...)```

Print the progress line for step ```i``` if one is due (see ```ProgressLog```). Returns \
whether a line was printed."""
function log_progress!(p::ProgressLog, i::Integer, values::Pair...)
  isdue(p, i) || return false
  print_line(p, i, values)
  return true
end

"""```needs_final(p::ProgressLog, i)```

Whether ```log_done!(p, i, …)``` would print, i.e. ```p``` is enabled and step ```i``` \
was not printed yet. Use it to skip computing an expensive value for the final line."""
needs_final(p::ProgressLog, i::Integer) = p.enabled && i != p.printed

"""```log_done!(p::ProgressLog, i, values...)```

Print the final line for step ```i```, unless that step was already printed or ```p``` \
is disabled."""
function log_done!(p::ProgressLog, i::Integer, values::Pair...)
  needs_final(p, i) || return nothing
  print_line(p, i, values)
  return nothing
end
