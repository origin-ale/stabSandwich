# Item 8 Verification: compare two `accuracy_dense.jl` outputs, sample by sample.
#
#   julia --project=. scripts/item8/accuracy_compare.jl old.jls new.jls
#
# For each sample: max |⟨O⟩ − dense| of each run, and the peak and summed layer-end max bond
# dimension of each. Then the totals: samples where new error > old error (and by how much),
# the largest error of each run, and the samples where the new layer-end bond dimension is
# above the old one at any layer end.
using Printf, Serialization

length(ARGS) == 2 || error("usage: accuracy_compare.jl old.jls new.jls")
old, new = deserialize(ARGS[1]), deserialize(ARGS[2])
old.N == new.N || error("N differs")
key(r) = (r.p, r.s)
oldby = Dict(key(r) => r for r in old.recs)
length(old.recs) == length(new.recs) && all(haskey(oldby, key(r)) for r in new.recs) ||
  error("the two files hold different samples")

@printf("%-6s %-3s %-6s %-11s %-11s %-9s %-9s %-9s %-9s\n", "p", "s", "doped",
        "err old", "err new", "peak old", "peak new", "Σbd old", "Σbd new")
worse = 0; worse_ratio = 1.0; bd_above = 0
for r in new.recs
  o = oldby[key(r)]
  o.ndoped == r.ndoped || error("doping differs at $(key(r))")
  eo, en = maximum(o.err), maximum(r.err)
  @printf("%-6.3f %-3d %-6d %-11.2e %-11.2e %-9d %-9d %-9d %-9d\n", r.p, r.s, r.ndoped, eo, en,
          maximum(o.bds), maximum(r.bds), sum(o.bds), sum(r.bds))
  if en > eo
    global worse += 1
    global worse_ratio = max(worse_ratio, en / max(eo, eps()))
  end
  any(r.bds .> o.bds) && (global bd_above += 1)
end
n = length(new.recs)
println()
@printf("N = %d, %d samples\n", new.N, n)
@printf("max |⟨O⟩ − dense|: old %.2e, new %.2e\n", maximum(maximum(r.err) for r in old.recs),
        maximum(maximum(r.err) for r in new.recs))
@printf("new error > old error in %d/%d samples (largest ratio %.1f)\n", worse, n, worse_ratio)
@printf("new layer-end max bond dim above old at some layer end: %d/%d samples\n", bd_above, n)
@printf("Σ over samples of peak bond dim: old %d, new %d\n", sum(maximum(r.bds) for r in old.recs),
        sum(maximum(r.bds) for r in new.recs))
