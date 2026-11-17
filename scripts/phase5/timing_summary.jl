# Phase 5: summarise per-sample times from `scripts/t11/compare_prob_tm.jl ... tsv=<file>` tables
# (old rerun vs new replay at the same thread count). Indicative only: the cluster is shared.
#
#   julia scripts/phase5/timing_summary.jl <table.tsv> [more.tsv ...]
#
# For each table, method and p: number of samples, mean old and new time, the median and
# minimum per-sample speed-up (old/new), and the slowest sample on each stack.

using Printf, Statistics

for path in ARGS
  lines = filter(!isempty, readlines(path))
  header = split(lines[1], '\t')
  col(name) = findfirst(==(name), header)
  im, ip, is, io, inew = col("method"), col("p"), col("sample"), col("time old (s)"), col("time new (s)")
  groups = Dict{Tuple{String, String}, Vector{Tuple{Int, Float64, Float64}}}()
  for l in lines[2:end]
    r = split(l, '\t')
    (r[io] == "-" || r[inew] == "-") && continue
    push!(get!(groups, (r[im], r[ip]), []), (parse(Int, r[is]), parse(Float64, r[io]), parse(Float64, r[inew])))
  end
  println("== ", path)
  @printf("%-8s %-5s %3s  %10s %10s  %8s %8s  %-22s %-22s\n", "method", "p", "n",
          "mean old", "mean new", "med ×", "min ×", "slowest old (sample)", "slowest new (sample)")
  for key in sort(collect(keys(groups)); by = k -> (k[1], parse(Float64, k[2])))
    v = groups[key]
    told, tnew = [x[2] for x in v], [x[3] for x in v]
    sp = told ./ tnew
    jo, jn = argmax(told), argmax(tnew)
    @printf("%-8s %-5s %3d  %10.3g %10.3g  %8.2f %8.2f  %-22s %-22s\n", key[1], key[2], length(v),
            mean(told), mean(tnew), median(sp), minimum(sp),
            @sprintf("%.3g s (%d)", told[jo], v[jo][1]), @sprintf("%.3g s (%d)", tnew[jn], v[jn][1]))
  end
  println()
end
