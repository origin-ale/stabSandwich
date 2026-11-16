# T11: per-sample old-vs-new comparison of the Phase 0 `Comptime_prob_tm.jl` runs and
# their replays (`replay_prob_tm.jl`). Run locally from the repository root:
#
#   julia --project=. scripts/t11/compare_prob_tm.jl old=<dir> new=<dir> [tsv=<file>]
#
#   old  directory with the Phase 0 output files, e.g.
#        output/runs/<Phase 0 id>/files/output (never output/ itself: the Phase 0
#        jobs share file names there)
#   new  directory with the replay's output files, e.g.
#        output/runs/<replay id>/files/output/t11/<Phase 0 id>
#   tsv  optional: also write the per-sample table to this file
#
# Compared for each method, p and sample of the replay (plan, "Old-vs-new regression
# design"; tolerances decided by the user 2026-11-10):
#   - inputs: initial bitstring and doped gates identical (a replay of other circuits
#     isn't a regression test);
#   - the run parameters in the two headers identical;
#   - CAMPS-PP: switch gate, k at the switch and PP stop gate identical; layers run in
#     CAMPS by both (the old bond-dimension block) |Δ⟨O⟩| ≤ 1e-10, later layers ≤ δ;
#   - CAMPS (`:camps`): every layer |Δ⟨O⟩| ≤ δ;
#   - the same layers recorded by both.
# δ = 1e-6 is never loosened here; a change needs the user's written approval.
# Resources (bond dimension, n. of Paulis) and times are printed side by side for the
# user to review, not pass/fail. Exit code 1 on any violation, 0 otherwise.

using Printf

const ATOL_PRESWITCH = 1e-10
const DELTA = 1e-6

args = Dict{String, String}()
for arg in ARGS
  kv = split(arg, "="; limit = 2)
  length(kv) == 2 && kv[1] in ("old", "new", "tsv") ||
    error("Arguments are old=<dir>, new=<dir> and optionally tsv=<file>; got `$arg`")
  args[kv[1]] = kv[2]
end
haskey(args, "old") && haskey(args, "new") || error("Need old=<dir> and new=<dir>")
olddir, newdir = args["old"], args["new"]

# == Parsers (layouts of initialize_output and the Utility.jl writers) =========
function read_header(file)
  line = first(eachline(file))
  startswith(line, "#") || error("$file has no header line")
  return Dict(String(kv[1]) => String(kv[2])
    for kv in (split(t, "="; limit = 2) for t in split(line[2:end])) if length(kv) == 2)
end

""" Sample ids in block order: `sample_ids` in a replay, 1:samples in a Phase 0 run. """
sample_ids(hdr) = haskey(hdr, "sample_ids") ? parse.(Int, split(hdr["sample_ids"], ",")) :
  collect(1:parse(Int, hdr["samples"]))

""" `save_full_samples` layout => (p, sample id) => layer => value. """
function read_full(file, ids)
  data = Dict{Tuple{Float64, Int}, Dict{Int, Float64}}()
  p = nothing; cur = nothing
  for line in eachline(file)
    if (m = match(r"^# p = (\S+), μ = ", line)) !== nothing
      p = parse(Float64, m[1]); cur = nothing
    elseif (m = match(r"^# sample (\d+)$", line)) !== nothing
      key = (p, ids[parse(Int, m[1])])
      haskey(data, key) && error("$file: $key listed twice")
      cur = data[key] = Dict{Int, Float64}()
    elseif !startswith(line, "#") && !isempty(strip(line))
      l, v = split(line, '\t')
      cur[parse(Int, l)] = parse(Float64, v)
    end
  end
  return data
end

""" `save_per_sample` layout => (p, sample id) => [columns...]. """
function read_per_sample(file, ids)
  data = Dict{Tuple{Float64, Int}, Vector{Float64}}()
  p = nothing
  for line in eachline(file)
    if (m = match(r"^# p = (\S+), μ = ", line)) !== nothing
      p = parse(Float64, m[1])
    elseif !startswith(line, "#") && !isempty(strip(line))
      t = split(line)
      data[(p, ids[parse(Int, t[1])])] = parse.(Float64, t[2:end])
    end
  end
  return data
end

""" `initstates_prob_tm.txt` => sample id => the rest of its line. """
read_initstates(file) = Dict(parse(Int, split(l)[1]) => join(split(l)[2:end], " ")
  for l in eachline(file) if !startswith(l, "#") && !isempty(strip(l)))

""" `save_doped_gates` layout => (p, sample id) => gate lines (block positions mapped to ids). """
function read_doped(file, ids)
  data = Dict{Tuple{Float64, Int}, Vector{String}}()
  p = nothing; cur = nothing
  for line in eachline(file)
    if (m = match(r"^# p = (\S+), μ = ", line)) !== nothing
      p = parse(Float64, m[1])
    elseif (m = match(r"^# sample (\d+): ", line)) !== nothing
      cur = data[(p, ids[parse(Int, m[1])])] = String[]
    elseif !startswith(line, "#") && !isempty(strip(line))
      push!(cur, strip(line))
    end
  end
  return data
end

# == Load ======================================================================
f(dir, name) = (path = joinpath(dir, name); isfile(path) || error("Missing $path"); path)

hdr_new = read_header(f(newdir, "initstates_prob_tm.txt"))
hdr_old = read_header(f(olddir, "initstates_prob_tm.txt"))
ids_new, ids_old = sample_ids(hdr_new), sample_ids(hdr_old)
methods = Symbol.(split(hdr_new["methods"], ","))
probs = parse.(Float64, split(hdr_new["magic_probs"], ","))
thl_exp = abs(round(Int, log10(parse(Float64, hdr_new["campspp_thl"]))))

violations = String[]
violation(msg) = (push!(violations, msg); println("VIOLATION: ", msg))

println("old: ", olddir)
println("new: ", newdir, " (replay_of = ", get(hdr_new, "replay_of", "?"),
  ", inject = ", get(hdr_new, "inject", "?"), ", stack = ", get(hdr_new, "stack", "?"), ")")

# Run parameters: everything except what a replay legitimately changes or adds
let skip = ("methods", "magic_probs", "samples", "sample_ids", "replay_of", "inject", "stack")
  for k in sort(collect(union(keys(hdr_old), keys(hdr_new))))
    k in skip && continue
    get(hdr_old, k, nothing) == get(hdr_new, k, nothing) ||
      violation("parameter $k: old $(get(hdr_old, k, "(none)")), new $(get(hdr_new, k, "(none)"))")
  end
end
issubset(methods, Symbol.(split(hdr_old["methods"], ","))) ||
  violation("replay methods $methods not all in the old run")
issubset(ids_new, ids_old) || violation("replay samples $ids_new not all in the old run")

# Inputs
init_old = read_initstates(f(olddir, "initstates_prob_tm.txt"))
init_new = read_initstates(f(newdir, "initstates_prob_tm.txt"))
doped_old = read_doped(f(olddir, "doped_prob_tm.txt"), ids_old)
doped_new = read_doped(f(newdir, "doped_prob_tm.txt"), ids_new)
for s in ids_new
  get(init_old, s, nothing) == get(init_new, s, nothing) || violation("sample $s: initial states differ")
  for p in probs
    get(doped_old, (p, s), nothing) == get(doped_new, (p, s), nothing) ||
      violation("p = $p, sample $s: doped gates differ")
  end
end

# == Compare ===================================================================
rows = Vector{Vector{String}}()
header = ["method", "p", "sample", "switch old", "switch new", "k old", "k new",
  "peak bd old", "peak bd new", "peak NP old", "peak NP new",
  "max|Δ| CAMPS layers", "max|Δ| other layers", "time old (s)", "time new (s)"]
fmt(x) = x === nothing ? "-" : x isa AbstractFloat ? @sprintf("%.3g", x) : string(x)
peak(d) = isempty(d) ? nothing : Int(maximum(values(d)))

for m in methods, p in probs
  evs_old = read_full(f(olddir, "evs_prob_tm_$(m)_$(thl_exp)_full.txt"), ids_old)
  evs_new = read_full(f(newdir, "evs_prob_tm_$(m)_$(thl_exp)_full.txt"), ids_new)
  bd_old = read_full(f(olddir, "resources_prob_tm_bd_full_$m.txt"), ids_old)
  bd_new = read_full(f(newdir, "resources_prob_tm_bd_full_$m.txt"), ids_new)
  t_old = read_per_sample(f(olddir, "comptime_full_prob_tm_$m.txt"), ids_old)
  t_new = read_per_sample(f(newdir, "comptime_full_prob_tm_$m.txt"), ids_new)
  if m == :campspp
    sw_old = read_per_sample(f(olddir, "campspp_switch_prob_tm.txt"), ids_old)
    sw_new = read_per_sample(f(newdir, "campspp_switch_prob_tm.txt"), ids_new)
    np_old = read_full(f(olddir, "resources_prob_tm_NP_full_campspp.txt"), ids_old)
    np_new = read_full(f(newdir, "resources_prob_tm_NP_full_campspp.txt"), ids_new)
  end

  for s in ids_new
    key = (p, s)
    tag = "$m p = $p, sample $s"
    if !(haskey(evs_old, key) && haskey(evs_new, key))
      violation("$tag: ⟨O⟩ missing (old $(haskey(evs_old, key)), new $(haskey(evs_new, key)))")
      continue
    end
    eo, en = evs_old[key], evs_new[key]
    keys(eo) == keys(en) ||
      violation("$tag: layers recorded differ (old $(sort(collect(keys(eo)))), new $(sort(collect(keys(en)))))")
    layers = sort(collect(intersect(keys(eo), keys(en))))

    # Layers run in CAMPS by both stacks, at 1e-10 (CAMPS-PP only, and only if the
    # switch matches; the old bond-dimension block lists exactly those layers)
    camps_layers = Int[]
    swo = swn = nothing
    if m == :campspp
      swo, swn = get(sw_old, key, nothing), get(sw_new, key, nothing)
      if swo === nothing || swn === nothing
        violation("$tag: switch record missing")
      elseif swo != swn
        violation("$tag: switch (gate, k, PP stop) old $(Int.(swo)), new $(Int.(swn))")
      else
        camps_layers = sort(collect(keys(bd_old[key])))
        keys(bd_old[key]) == keys(bd_new[key]) ||
          violation("$tag: CAMPS layers differ despite the same switch gate")
      end
    end

    dpre = dpost = nothing
    for l in layers
      d = abs(eo[l] - en[l])
      if l in camps_layers
        dpre = max(something(dpre, 0.0), d)
        d > ATOL_PRESWITCH && violation(@sprintf("%s, layer %d (CAMPS): |Δ⟨O⟩| = %.3e > %.0e", tag, l, d, ATOL_PRESWITCH))
      else
        dpost = max(something(dpost, 0.0), d)
        d > DELTA && violation(@sprintf("%s, layer %d: |Δ⟨O⟩| = %.3e > δ = %.0e", tag, l, d, DELTA))
      end
    end

    push!(rows, fmt.([m, p, s,
      swo === nothing ? nothing : Int(swo[1]), swn === nothing ? nothing : Int(swn[1]),
      swo === nothing ? nothing : Int(swo[2]), swn === nothing ? nothing : Int(swn[2]),
      peak(get(bd_old, key, Dict())), peak(get(bd_new, key, Dict())),
      m == :campspp ? peak(get(np_old, key, Dict())) : nothing,
      m == :campspp ? peak(get(np_new, key, Dict())) : nothing,
      dpre, dpost,
      haskey(t_old, key) ? t_old[key][1] : nothing, haskey(t_new, key) ? t_new[key][1] : nothing]))
  end
end

# == Report ====================================================================
widths = [maximum(length, [header[j]; [r[j] for r in rows]]) for j in eachindex(header)]
println()
println(join(rpad.(header, widths), "  "))
for r in rows
  println(join(rpad.(r, widths), "  "))
end
if haskey(args, "tsv")
  open(args["tsv"], "w") do io
    println(io, join(header, '\t'))
    foreach(r -> println(io, join(r, '\t')), rows)
  end
end

println()
println("Compared $(length(rows)) (method, p, sample) cases. Tolerances: CAMPS layers of CAMPS-PP ",
  ATOL_PRESWITCH, ", everything else δ = ", DELTA, ".")
if isempty(violations)
  println("PASS: no violations.")
else
  println("FAIL: $(length(violations)) violation(s).")
  exit(1)
end
