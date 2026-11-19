# Item 8, decision 1: analyse `tau_measure.jl` output. Run from the repository root:
#
#   julia --project=. scripts/item8/tau_analyse.jl out.jls [out2.jls ...]
#
# Per file: a decade histogram of σ/σ_max (identity candidate, the 19 gate candidates, and
# spectra taken after the `truncated` flag), the widest empty gap among values in [1e-18, 1),
# and for each τ = 1e-16 … 1e-4 the number of values within a factor of 100 of τ (item 8's
# Verification criterion: none may be) and the number of bond searches in which some gate
# lowers the τ-rank below the identity's (where `:chi3` with that τ would act). Last, the
# stored bonds whose dimension exceeds their numerical rank at τ = 1e-10 (zero padding).
using Serialization, Printf
# Files written before the switch to plain Vector{Float64} hold NDTensors vectors
try
  @eval using ITensors
catch
end

const TAU = 1e-10

for f in ARGS
  d = deserialize(f)
  recs = [merge(r, (rel = Float64[x for x in r.rel],)) for r in d.recs]
  println("\n==== $(basename(f)): N=$(d.N), $(length(recs)) candidate spectra, ",
          length(unique((r.p, r.s, r.layer, r.bond) for r in recs)), " (run, layer, bond) searches")
  for (lab, sel) in (("identity (two-site SVD of the stored bond)", r -> r.cand == 0),
                     ("19 gate candidates", r -> r.cand != 0),
                     ("after the truncated flag", r -> r.truncated))
    vals = reduce(vcat, [r.rel for r in recs if sel(r)]; init = Float64[])
    isempty(vals) && (println("  $lab: none"); continue)
    println("  $lab: $(length(vals)) values, exact zeros $(count(==(0.0), vals))")
    c = count(x -> 0 < x <= 1e-18, vals)
    c > 0 && @printf("    (0, 1e-18]  %d\n", c)
    for e in -17:0
      c = count(x -> 10.0^(e - 1) < x <= 10.0^e, vals)
      c > 0 && @printf("    (1e%d, 1e%d]  %d\n", e - 1, e, c)
    end
  end
  lv = sort(unique(log10.(filter(x -> 1e-18 <= x < 1, reduce(vcat, [r.rel for r in recs])))))
  if length(lv) > 1
    gaps = diff(lv); k = argmax(gaps)
    @printf("  widest empty gap in [1e-18, 1): (%.3g, %.3g), %.1f decades\n",
            10^lv[k], 10^lv[k+1], gaps[k])
  end

  searches = Dict{Any, Vector{Any}}()
  for r in recs
    push!(get!(searches, (r.p, r.s, r.layer, r.bond), Any[]), r)
  end
  println("  τ       vals in [τ/100,100τ]  searches where a gate lowers the τ-rank")
  for e in -16:-4
    τ = 10.0^e
    near = sum(r -> count(x -> τ / 100 <= x <= 100τ, r.rel), recs)
    act = count(values(searches)) do rs
      rk = Dict(r.cand => count(>(τ), r.rel) for r in rs)
      minimum(rk[c] for c in keys(rk) if c != 0) < rk[0]
    end
    @printf("  1e%-4d  %8d              %d / %d\n", e, near, act, length(searches))
  end

  ids = [r for r in recs if r.cand == 0]
  pad = [r for r in ids if count(>(TAU), r.rel) < r.dim]
  println("  stored bonds whose dimension exceeds their rank at τ = $TAU: ",
          length(pad), " / ", length(ids))
  for r in pad[1:min(end, 6)]
    println("    p=$(r.p) s=$(r.s) layer=$(r.layer) bond=$(r.bond) dim=$(r.dim) rel σ = ",
            round.(r.rel; sigdigits = 3))
  end
  println("  runs:")
  foreach(s -> println("    ", s), d.summary)
end
