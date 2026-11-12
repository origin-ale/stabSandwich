# Random Clifford-angle rotation circuits at N = 46: the CampsPP dynamics against the
# stabilizer reference (T6, large N). Expected values: StabilizerReference.jl, which
# implements CONVENTIONS.md §2 and is checked against the dense reference in
# CampsPP/test/DisentangleCAMPS/TestCliffordLargeN.jl.
#
# Unlike the undoped XXZ test (a bit permutation), these circuits spread the stabilizers
# over all 46 qubits, so the CAMPS frame is a dense Clifford and PP's Pauli strings are high
# weight with UInt96 keys. The hand-over tests evaluate a 46-qubit Pauli sum on such a CAMPS
# through the old pp → cmps converter and `expectation`.

isdefined(@__MODULE__, :StabilizerReference) ||
  include(joinpath(@__DIR__, "DisentangleCAMPS", "StabilizerReference.jl"))
const SR = StabilizerReference

"String over IXYZ → Hermitian PauliOperator (phase byte 0)."
strop(str) = PauliOperator(0x0, [c in "XY" for c in str], [c in "ZY" for c in str])

"⟨O⟩ of a `pp.PauliSum` on the stabilizer state with generators `gens` (term bits by §3)."
function stab_expect(gens, O::pp.PauliSum)
  S = Stabilizer(gens)
  return sum(coeff * SR.stab_expect(S, strop(termstring(term, O.nqubits))) for (term, coeff) in O.terms)
end

"A signed PauliOperator as a `pp.PauliString` with coefficient `c` × its sign."
function pp_term(g::PauliOperator, c)
  N = length(g)
  qs = [j for j in 1:N if xbit(g)[j] || zbit(g)[j]]
  syms = [xbit(g)[j] ? (zbit(g)[j] ? :Y : :X) : :Z for j in qs]
  sign = g.phase[] == 0x2 ? -1 : 1
  return pp.PauliString(N, syms, qs, sign * c)
end

"Reference ⟨O⟩ at t = 0 and after each gate index in `record_at`; also the generators there."
function stab_dynamics(N, ones, gates, phases, O, record_at)
  gens = SR.basis_generators(N, ones)
  evs = [stab_expect(gens, O)]
  snaps = []
  for (i, (g, ϕ)) in enumerate(zip(gates, phases))
    SR.rotate_generators!(gens, strop(gatestring(g, N)), ϕ)
    if i in record_at
      push!(evs, stab_expect(gens, O))
      push!(snaps, copy(gens))
    end
  end
  return evs, snaps
end

@testset "random Clifford circuits at N = 46 vs stabilizer reference" begin
  N = 46
  slow = get(ENV, "STABSANDWICH_SLOW", "") != ""
  M = slow ? 200 : 60
  le = collect(M÷4:M÷4:M)
  _, ones = domainwallstate(Xoshiro(100), N, 0.6)
  rng = Xoshiro(1301)
  # Angles: the Clifford branches that agree with exp(iϕP) in the old code (CONVENTIONS.md §2).
  ϕs = rand(rng, [π/4, π/4, -3π/4, 0.0, π], M)
  ϕs[1] = π/4                                # so a hand-over before gate 1 is visible
  gates, phases = rotation_circuit(rng, ϕs, N)
  @test count(abs.(phases) .≈ π/4) >= M / 4

  # Hand-over points. 0 and π gates are the identity, so a hand-over one gate off is
  # invisible next to them: use points s0 where gates s0 and s0+1 are both nontrivial.
  nontrivial(i) = !(1 <= i <= M) || !any(isapprox(abs(phases[i]), c; atol = 1e-12) for c in (0.0, π))
  visible(s0) = nontrivial(s0) && nontrivial(s0 + 1)
  up(c) = findfirst(visible, c:M) + c - 1
  down(c) = c - findfirst(visible, c:-1:0) + 1
  handovers = unique([0, up(7), le[1], up(le[1]), up(le[2] + 1), down(M - 1), M])
  # (le[1] is kept as is, to test a hand-over exactly at a record point; M leaves nothing for PP.)
  @test all(visible, filter(s0 -> s0 ∉ (le[1], M), handovers))
  @test length(handovers) >= 6

  # Observable: transferred magnetization plus, for each record point, two of the reference
  # generators there (⟨·⟩ = 1 at that point, in general 0 elsewhere), with distinct weights.
  # Built from a reference run first, so every record point has a nontrivial value.
  O0 = transferredmagnetization(N, ones)
  _, snaps = stab_dynamics(N, ones, gates, phases, O0, le)
  O = deepcopy(O0)
  for (l, gens) in enumerate(snaps), j in 1:2
    O += pp_term(gens[7j], 0.5^l * j)
  end
  ref, _ = stab_dynamics(N, ones, gates, phases, O, le)
  @test length(ref) == length(le) + 1
  @test all(abs(ref[l+1] - stab_expect(snaps[l], O0)) > 0.1 for l in eachindex(le))
  # The reference is sensitive to the sign of the angles and to the gate order.
  @test maxabsdiff(ref, first(stab_dynamics(N, ones, gates, -phases, O, le))) > 0.1
  @test maxabsdiff(ref, first(stab_dynamics(N, ones, reverse(gates), reverse(phases), O, le))) > 0.1
  # The generators are spread over the register, not a bit permutation.
  @test sum(count(xbit(g) .| zbit(g)) for g in snaps[end]) / N > N / 3

  ψ0 = computationalcamps(N, ones)

  @testset "camps_circuit_dynamics" begin
    χ_camps = 64
    ψ, k, s, evs, bds = camps_circuit_dynamics(ψ0, gates, phases, χ_camps, O; layer_ends = le, track = true)
    @test s == M
    @test maxabsdiff(evs, ref) <= EXACT_ATOL
    @test k == 0 && all(==(1), bds)            # frame-dependent: every gate went to the frame
  end

  @testset "pauliprop_circuit_dynamics (PP alone)" begin
    thl_pp, Nmax_pp = 1e-10, 10^6
    s, evs, nterms = pauliprop_circuit_dynamics(ones, gates, phases, thl_pp, Nmax_pp, O;
      layer_ends = le, track = true)
    @test s == M
    @test maxabsdiff(evs, ref) <= EXACT_ATOL
    # A Clifford rotation maps each string to one string; the cos(π/2) ≈ 6e-17 branches
    # fall below the threshold, so the sum keeps its size.
    @test all(==(length(O)), nterms)
  end

  @testset "PP continues from a CAMPS state after $s0 gates" for s0 in handovers
    thl_pp, Nmax_pp = 1e-10, 10^6
    ψ, k = evolve(computationalcamps(N, ones), s0, gates, phases)
    @test k == 0                               # frame-dependent
    s, evs = pauliprop_circuit_dynamics(ψ, s0, gates, phases, thl_pp, Nmax_pp, O; layer_ends = le)
    later = [l for l in eachindex(le) if le[l] > s0]
    @test s == M
    @test length(evs) == length(later)
    isempty(later) || @test maxabsdiff(evs, ref[later .+ 1]) <= EXACT_ATOL
  end

  @testset "campspp_circuit_dynamics, χ = $χ_campspp" for χ_campspp in (1, 64)
    thl_campspp, Nmax_campspp = 1e-10, 10^6
    evs, s, _, _, _, _, s_switch, k_switch = campspp_circuit_dynamics(
      ψ0, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, O; layer_ends = le, track = true)
    @test s == M
    @test maxabsdiff(evs, ref) <= EXACT_ATOL
    # Frame-dependent: the bond dimension is 1 throughout, so χ = 1 hands over before the
    # first gate (the initial state already has bond dimension ≥ χ) and χ = 64 never does.
    # A Clifford-only circuit cannot force a switch in mid-circuit, so the glue between the
    # two stages (which state and which gate PP starts from) is not tested here; that is
    # covered at small N in TestTimeEvolution.jl (χ = 2, 4), and the hand-over itself on a
    # dense N = 46 frame by the "PP continues from a CAMPS state" tests above.
    @test s_switch == (χ_campspp == 1 ? 0 : M) && k_switch == 0
  end

  @testset "campssrc_circuit_dynamics ($crit, $strat)" for (crit, strat) in
      [(:chi3, :snake), (:entangle, :brickwork)]
    thl_camps = 1e-12
    ψ, s, evs = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
      criterion = crit, strategy = strat, layer_ends = le)
    @test s == M
    @test maxabsdiff(evs, ref) <= EXACT_ATOL
  end
end
