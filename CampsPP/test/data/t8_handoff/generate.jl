# T8 reference data (plan, "Converter at large sums"): real CAMPS-PP hand-off sums at N = 46
# and the OLD code's ⟨O⟩ on a real hand-off CAMPS state, saved for the new code to reproduce.
#
# One sample of the production ensemble (Comptime_prob_tm.jl construction, rng_seed = 100)
# at p = P_DOPE. The state is the CAMPS after its switch gate s (bond dimension ≥ 64). Each
# sum is the transferred magnetization propagated by PP (thl = 1e-10) through gates s+1…e,
# i.e. exactly what CAMPS-PP hands to the converter when it records at gate e; e is chosen
# so the sum has about the target number of terms.
#
# Usage: julia --project=. t8_handoff.jl <outdir> [p] [sample] [targets…]

using CampsPP, CampsPP.DisentangleCAMPS
using ITensors, ITensorMPS
import PauliPropagation as pp
import CliffordMPS as cmps
using QuantumClifford
using Random: MersenneTwister
using Printf, LinearAlgebra, Strided
Strided.disable_threads()
BLAS.set_num_threads(Threads.nthreads())

outdir = ARGS[1]
const P_DOPE = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 0.06
const SAMPLE = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1
targets = length(ARGS) >= 4 ? parse.(Int, ARGS[4:end]) : [100, 500, 1000]
const N, μ, rng_seed, χ, thl = 46, 0.6, 100, 64, 1e-10
const magic_phase = 3π/16
mkpath(outdir)

# Spec bit rule (CONVENTIONS.md §3), independent of the converter under test.
termstring(term, N) = String(["IXYZ"[Int((term >> (2 * (j - 1))) & 0x3) + 1] for j in 1:N])
function fromstrings(N, rows)
  S = pp.PauliSum(N)
  for (str, c) in rows
    qs = [j for j in 1:N if str[j] != 'I']
    pp.add!(S, pp.PauliString(N, [Symbol(str[j]) for j in qs], qs, c))
  end
  return S
end
strop(str) = PauliOperator(0x0, [c in "XY" for c in str], [c in "ZY" for c in str])

_, ones = domainwallstate(MersenneTwister(rng_seed + SAMPLE), N, μ)
gates, phases_ref = xxz_circuit(π/4, π/4, N/2, N)
mask = doping_mask(MersenneTwister(hash((rng_seed + SAMPLE, P_DOPE))), phases_ref, P_DOPE)
phases = xy_magic(phases_ref, mask; magicphase = magic_phase)
doped = findall(phases .!= phases_ref)
M = length(gates)
O = transferredmagnetization(N, ones)

# CAMPS up to the switch (the same loop CAMPS-PP runs; untruncated).
ψ, k, s, _ = camps_circuit_dynamics(computationalcamps(N, ones), gates, phases, χ, O)
bd = DisentangleCAMPS.bonddim(ψ)
@info "switch" s k bd
s < M || error("sample does not switch")

const PPALONE = get(ENV, "MODE", "") == "ppalone"
# Default: the real hand-off (gates s+1…e). ppalone: PP alone from t = 0 (gates 1…e), for sums
# larger than this sample's hand-off reaches; evaluated on the same state.
first_gate = PPALONE ? 1 : s + 1
handoff(e) = pp.propagate(gates[first_gate:e], O, -2 .* phases[first_gate:e]; min_abs_coeff = thl)

isfile(joinpath(outdir, "state.txt")) || open(joinpath(outdir, "state.txt"), "w") do io
  println(io, "# T8 hand-off state (old code, CliffordMPS frame). Rebuild: basis state with the")
  println(io, "# excited qubits below; xxz_circuit(π/4, π/4, 23, 46); replace the XX and YY angles of")
  println(io, "# the doped gates below by 3π/16; apply gates 1…s untruncated. ⟨O⟩ is frame-invariant.")
  println(io, "N\t$N")
  println(io, "p\t$P_DOPE")
  println(io, "sample\t$SAMPLE (Comptime_prob_tm.jl ensemble, rng_seed = $rng_seed, μ = $μ)")
  println(io, "s\t$s")
  println(io, "k\t$k (old code, frame-dependent)")
  println(io, "bonddim\t$bd (old code, frame-dependent)")
  println(io, "excited\t", join(ones, ' '))
  println(io, "doped_gates\t", join(doped, ' '))
end

summary = []
for target in targets
  # Grow the end gate one layer at a time (layer l ends at gate 3(N−1)·l) and stop at the
  # first layer end whose sum reaches the target; then bisect inside that layer. This never
  # propagates far past the target size (PP alone through the whole circuit can reach
  # 10⁷–10⁸ terms at p = 0.06).
  lstep = 3 * (N - 1)
  ends = [e for e in lstep:lstep:M if e >= first_gate]
  prev, hi = first_gate - 1, 0
  # Explicit end gate (ENV["E"]): skip the search. Sizes are not monotone in e inside an
  # XX-YY-ZZ triple (after XX the sum can be 50× larger than after the full triple), so a
  # bisection can land on a huge mid-triple sum.
  if haskey(ENV, "E")
    ends = Int[]; prev = hi = lo0 = parse(Int, ENV["E"])
  end
  for e in ends
    n = length(handoff(e))
    println("  target $target: gate $e → $n terms"); flush(stdout)
    n >= target && (hi = e; break)
    prev = e
  end
  hi == 0 && (@warn "target $target not reached by the end of the circuit"; continue)
  lo = haskey(ENV, "E") ? hi : max(prev + 1, first_gate)
  while lo < hi
    mid = (lo + hi) ÷ 2
    length(handoff(mid)) >= target ? (hi = mid) : (lo = mid + 1)
  end
  e = lo
  S = handoff(e)
  rows = sort([(termstring(t, N), c) for (t, c) in S.terms])
  @assert all(length(r[1]) == N for r in rows)
  # Serialisation round trip: the file rebuilds the same sum.
  @assert fromstrings(N, rows).terms == S.terms

  GC.gc()
  t_mpo = @elapsed ev_mpo = cmps.expectation(ψ, S)                  # old converter + MPO route
  t_terms = @elapsed ev_terms = sum(c * cmps.expectation(ψ, strop(str)) for (str, c) in rows)
  @info "sum" target e length(S) ev_mpo ev_terms t_mpo t_terms
  path = joinpath(outdir, @sprintf("sum_%04d.txt", length(S)))
  open(path, "w") do io
    println(io, "# T8 hand-off sum: transferred magnetization propagated by PP (thl = $thl) through")
    println(io, "# gates $(first_gate)…e of the circuit in state.txt", PPALONE ? " (PP alone from t = 0; not the physical hand-off for this state)." : " (s+1…e: the real hand-off).")
    println(io, "# Rows: Pauli string (character j = qubit j), coefficient. ⟨O⟩ is on the state in state.txt.")
    println(io, "# e\t$e")
    println(io, "# nterms\t$(length(S))")
    println(io, "# old_expectation_mpo\t", repr(ev_mpo))
    println(io, "# old_expectation_termwise\t", repr(ev_terms))
    @printf(io, "# old_time_mpo_s\t%.3f (local, single @elapsed, includes compilation for this size)\n", t_mpo)
    for (str, c) in rows
      println(io, str, '\t', repr(c))
    end
  end
  push!(summary, (; target, e, n = length(S), ev_mpo, ev_terms, t_mpo, t_terms))
end
println("p=$P_DOPE sample=$SAMPLE s=$s k=$k bd=$bd")
for r in summary
  @printf("target %5d: e=%4d nterms=%5d  ⟨O⟩mpo=%.15g  termwise Δ=%.1e  t_mpo=%.1fs t_terms=%.1fs\n",
    r.target, r.e, r.n, real(r.ev_mpo), abs(r.ev_mpo - r.ev_terms), r.t_mpo, r.t_terms)
end
