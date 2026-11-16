# T11 (optional): old vs new `Randrot_complexity` curves from `randrot_complexity.jl`.
#
#   julia scripts/t11/compare_randrot.jl <old.dat> <new.dat>
#
# M₂/N (frame-invariant): per point z = (new − old) / √(err_old² + err_new²). The plan asks
# for agreement "within the error bars"; here a point fails if |z| > 3. With 11 independent
# points, a 2σ cut would flag at least one point by chance about 40% of the time. Points
# with 2 < |z| ≤ 3 are listed for review.
# EE/N (frame-dependent): printed side by side, not pass/fail.
# Exit 1 if any M₂ point fails.

using Printf

function read_dat(file)
  rows = [parse.(Float64, split(l)) for l in eachline(file) if !startswith(l, "#") && !isempty(strip(l))]
  return reduce(hcat, rows)'   # columns: t/N, EE, err EE, M2, err M2, time
end

length(ARGS) == 2 || error("usage: compare_randrot.jl <old.dat> <new.dat>")
old, new = read_dat(ARGS[1]), read_dat(ARGS[2])
size(old) == size(new) && old[:, 1] == new[:, 1] || error("The two files have different t/N points")

println(first(eachline(ARGS[1])))
println(first(eachline(ARGS[2])))
@printf("%5s | %8s %8s %8s %8s %6s | %8s %8s %8s %8s | %6s %6s\n", "t/N",
  "M2 old", "± old", "M2 new", "± new", "z", "EE old", "± old", "EE new", "± new", "s old", "s new")
fails = 0
review = 0
for i in axes(old, 1)
  z = (new[i, 4] - old[i, 4]) / sqrt(old[i, 5]^2 + new[i, 5]^2)
  flag = abs(z) > 3 ? "  FAIL" : abs(z) > 2 ? "  review" : ""
  global fails += abs(z) > 3
  global review += 2 < abs(z) <= 3
  @printf("%5.1f | %8.4f %8.4f %8.4f %8.4f %6.2f | %8.4f %8.4f %8.4f %8.4f | %6.1f %6.1f%s\n",
    old[i, 1], old[i, 4], old[i, 5], new[i, 4], new[i, 5], z,
    old[i, 2], old[i, 3], new[i, 2], new[i, 3], old[i, 6], new[i, 6], flag)
end
println()
println("M₂/N: $(size(old, 1)) points, $fails with |z| > 3, $review with 2 < |z| ≤ 3.")
println(fails == 0 ? "PASS" : "FAIL")
exit(fails == 0 ? 0 : 1)
