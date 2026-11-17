# Phase 5: generate real N = 46 Pauli sums of about 10³, 10⁴ and 10⁵ terms for the converter and
# `expectation` benchmarks. PP alone from t = 0 (gates 1…e of the T8 sample's doped circuit,
# thl = 1e-10), as for T8's `sum_0538.txt`: real PP term structure, but ⟨O⟩ on the T8 state has
# no physical meaning (benchmark input only). Run from the repository root, under a timeout:
#
#   julia --project=. scripts/phase5/gen_sums.jl [E=1620,1755] [scan=2298:3:2427] [targets=1000,10000,100000] [maxterms=200000] [thl=1e-10]
#
# With E: one sum per listed end gate. With scan=lo:step:hi: print the term count at each gate in
# the range, save nothing, stop past `maxterms` (to choose E). Without: scan the layer ends after layer 11 (layer l ends at
# gate 135·l), save the first sum reaching each target, and stop once a sum exceeds `maxterms`.
# Layer ends only: inside an XX-YY-ZZ triple the term count is not monotone (T8 README).
# thl: PP coefficient threshold (default 1e-10, the production value); a larger one gives
# smaller real PP outputs at the same gate. First checks that e = 1504 reproduces T8's sum_0538.txt exactly. Output: output/phase5/sums/.

using CampsPP
include("common.jl")

args = Dict(split(a, "="; limit = 2) |> kv -> kv[1] => kv[2] for a in ARGS)
for k in keys(args)
  k in ("E", "scan", "targets", "maxterms", "thl") || error("Unrecognised argument `$k` (expected E, scan, targets, maxterms, thl)")
end
targets = parse.(Int, split(get(args, "targets", "1000,10000,100000"), ','))
maxterms = parse(Int, get(args, "maxterms", "200000"))
thl = parse(Float64, get(args, "thl", string(THL)))

st = read_state()
N = st.N
gates, phases = sample_circuit(st, xxz_circuit)
M = length(gates)
O = transferredmagnetization(N, st.excited)
ppalone(e; thl = thl) = pp.propagate(gates[1:e], O, -2 .* phases[1:e]; min_abs_coeff = thl)

# Sanity check against the T8 data: same construction, same sum.
ref538, _, _ = read_sum(joinpath(T8DIR, "sum_0538.txt"))
t1504 = @elapsed S1504 = ppalone(1504; thl = THL)
S1504.terms == ref538.terms || error("e = 1504 does not reproduce sum_0538.txt")
flushprintln("check: e = 1504 reproduces sum_0538.txt ($(length(S1504)) terms, $(round(t1504; digits = 1)) s incl. compilation)")

function save(e, S)
  mkpath(SUMDIR)
  rows = sort([(termstring(t, N), c) for (t, c) in S.terms])
  @assert ppsum_from_rows(N, rows).terms == S.terms   # file round trip
  path = joinpath(SUMDIR, @sprintf("sum_%06d.txt", length(S)))
  open(path, "w") do io
    println(io, "# Phase 5 benchmark sum: transferred magnetization propagated by PP alone (thl = $thl)")
    println(io, "# through gates 1…e of the T8 sample's circuit (CampsPP/test/data/t8_handoff/state.txt).")
    println(io, "# Not a physical hand-off: benchmark input only. Rows: Pauli string (character j = qubit j), coefficient.")
    println(io, "# e\t$e")
    println(io, "# thl\t$thl")
    println(io, "# nterms\t$(length(S))")
    for (str, c) in rows
      println(io, str, '\t', repr(c))
    end
  end
  flushprintln("  saved $path")
end

if haskey(args, "scan")
  lo, stp, hi = parse.(Int, split(args["scan"], ':'))
  for e in lo:stp:hi
    1 <= e <= M || error("scan gate $e outside 1…$M")
    local t = @elapsed n = length(ppalone(e))
    flushprintln(@sprintf("gate %4d → %7d terms (%.1f s)", e, n, t))
    n > maxterms && break
  end
elseif haskey(args, "E")
  for e in parse.(Int, split(args["E"], ','))
    1 <= e <= M || error("E = $e outside 1…$M")
    local t = @elapsed S = ppalone(e)
    flushprintln(@sprintf("gate %4d → %7d terms (%.1f s)", e, length(S), t))
    save(e, S)
  end
else
  pending = sort(targets)
  for e in 135 * 12:135:M
    local t = @elapsed S = ppalone(e)
    n = length(S)
    flushprintln(@sprintf("gate %4d → %7d terms (%.1f s)", e, n, t))
    reached = filter(<=(n), pending)
    isempty(reached) || save(e, S)
    filter!(>(n), pending)
    (isempty(pending) || n > maxterms) && break
  end
  isempty(pending) || @warn "targets not reached: $pending"
end
