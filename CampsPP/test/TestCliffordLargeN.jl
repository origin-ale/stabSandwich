# Random Clifford-angle rotation circuits at N = 46: the CampsPP dynamics against the
# stabilizer reference (T6, large N). Expected values: StabilizerReference.jl, which
# implements CONVENTIONS.md §2 and is checked against the dense reference in
# CampsPP/test/DisentangleCAMPS/TestCliffordLargeN.jl.
#
# Unlike the undoped XXZ test (a bit permutation), these circuits spread the stabilizers
# over all 46 qubits, so the CAMPS frame is a dense Clifford and PP's Pauli strings are high
# weight with keys wider than 64 bits (UInt128 on PP ≥ 0.9, UInt96 on 0.7.3). The hand-over tests evaluate a 46-qubit Pauli sum on such a CAMPS
# through the pp → QubitPauliSum converter and the MagicTensors `expectation` seam (T8).

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
  # Angles: the Clifford branches the old code got right (CONVENTIONS.md §2), kept from T-pre;
  # the other special angles are covered at N = 46 in DisentangleCAMPS/TestCliffordLargeN.jl.
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

  @testset "camps_circuit_dynamics, $(stackname(T))" for T in STATE_TYPES
    ψ0 = computationalcamps(T, N, ones)
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

  @testset "PP continues from a CAMPS state after $s0 gates, $(stackname(T))" for T in STATE_TYPES,
      s0 in handovers
    thl_pp, Nmax_pp = 1e-10, 10^6
    ψ, k = evolve(computationalcamps(T, N, ones), s0, gates, phases)
    @test k == 0                               # frame-dependent
    s, evs = pauliprop_circuit_dynamics(ψ, s0, gates, phases, thl_pp, Nmax_pp, O; layer_ends = le)
    later = [l for l in eachindex(le) if le[l] > s0]
    @test s == M
    @test length(evs) == length(later)
    isempty(later) || @test maxabsdiff(evs, ref[later .+ 1]) <= EXACT_ATOL
  end

  @testset "campspp_circuit_dynamics, χ = $χ_campspp, $(stackname(T))" for T in STATE_TYPES,
      χ_campspp in (1, 64)
    ψ0 = computationalcamps(T, N, ones)
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

  @testset "campssrc_circuit_dynamics ($crit, $strat), $(stackname(T))" for T in STATE_TYPES,
      (crit, strat) in [(:chi3, :snake), (:entangle, :brickwork)]
    ψ0 = computationalcamps(T, N, ones)
    thl_camps = 1e-12
    ψ, s, evs = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
      criterion = crit, strategy = strat, layer_ends = le)
    @test s == M
    @test maxabsdiff(evs, ref) <= EXACT_ATOL
  end
end

# -- Kicked Ising at its Clifford point, N = 46 (CONVENTIONS.md §2 "Kicked Ising") ----------
# ε = 0: every gate is Clifford. Reference: the spec circuit (KimSpec.jl) on the stabilizer
# reference, from the Néel state. ⟨Z_j⟩ itself is ±1 at t = 0 and 0 after layers 1 … 4 (the
# zero window of TestCircuits.jl; at ε = 0 and N = 46 it lasts beyond layer 30), so a second
# observable adds, for each layer l, two reference generators of the state after layer l
# (⟨·⟩ = 1 there, 0 at the other layers here) with weights 0.5^l and 0.5^l/3: its values are
# not 0, and they depend on the circuit. It is built from the spec (Z_j as a PP string), not
# from `local_z`, so the reference checks run before the KIM code exists. (A flipped ZZ sign
# is invisible to it, as to any physical check here; the structure tests in TestCircuits.jl
# catch that.) This is also the value check of the Clifford point for the KIM drivers: at
# N = 6 and ε = 0, ⟨Z_j⟩ is 0 at every layer the driver tests record (TestTimeEvolution.jl).
@testset "KIM at ε = 0, N = 46 vs stabilizer reference" begin
  N, t = 46, 4
  L = kim_spec_L(N)
  M = L * t
  le = [L * l for l in 1:t]
  ones = kim_neel(N)

  # Reference generators after each layer of the spec circuit (`wrap = false`: open chain)
  function kim_snaps(; wrap = true)
    # (`local`: names of the enclosing scope would otherwise be shared by this closure)
    local strings, ϕs = kim_spec_circuit(0.0, t, N; wrap)
    local Lw = length(strings) ÷ t
    local gens = SR.basis_generators(N, ones)
    local out = []
    for l in 1:t
      local r = (l - 1) * Lw + 1:l * Lw
      gens = SR.run_generators(gens, strop.(strings[r]), ϕs[r])
      push!(out, gens)
    end
    return out
  end
  snaps = kim_snaps()
  open_snaps = kim_snaps(wrap = false)
  ref_of(O, gs) = [stab_expect(SR.basis_generators(N, ones), O); [stab_expect(g, O) for g in gs]]
  spec_z(j) = pp.PauliSum(N) + pp.PauliString(N, [:Z], [j], 1.0)

  sites = (1, 2, 23, 46)
  observables = Dict(map(sites) do j
    Oaug = spec_z(j)
    for l in 1:t
      # a generator from around the chain, and generator 1, whose qubit sits on the wrap
      # bond (with t = 4 layers, generators far from (1, N) cannot see it)
      Oaug += pp_term(snaps[l][mod1(j + 7l, N)], 0.5^l)
      Oaug += pp_term(snaps[l][1], 0.5^l / 3)
    end
    j => (Oaug, ref_of(spec_z(j), snaps), ref_of(Oaug, snaps))
  end)

  @testset "references, j = $j" for j in sites
    Oaug, refz, refaug = observables[j]
    @test refz == [j in ones ? -1.0 : 1.0; zeros(t)]          # the zero window
    @test all(abs.(refaug[2:end]) .> 0.05)
    # Negative: the open chain (wrap bond dropped) changes the values (by 0.5²/3 at layer 2)
    @test maxabsdiff(refaug, ref_of(Oaug, open_snaps)) > 0.05
  end

  gates, phases = kicked_ising_circuit(0.0, t, N)
  @test length(gates) == M
  for j in sites
    Oaug, refz, refaug = observables[j]
    Oz = local_z(N, j)

    @testset "pauliprop_circuit_dynamics, j = $j" begin
      thl_pp, Nmax_pp = 1e-10, 10^6
      for (O, ref) in ((Oz, refz), (Oaug, refaug))
        s, evs, nterms = pauliprop_circuit_dynamics(ones, gates, phases, thl_pp, Nmax_pp, O;
          layer_ends = le, track = true)
        @test s == M
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
        # Each Clifford rotation maps a string to one string, distinct strings to distinct
        # strings; the cos(π/2) ≈ 6e-17 branches fall below thl_pp. So local_z stays 1
        # term, and the sum keeps its size.
        @test all(==(length(O)), nterms)
      end
    end

    @testset "campspp_circuit_dynamics, j = $j, $(stackname(T))" for T in STATE_TYPES
      χ_campspp, thl_campspp, Nmax_campspp = 64, 1e-10, 10^6
      for (O, ref) in ((Oz, refz), (Oaug, refaug))
        evs, s, bds, _, _, _, s_switch, k_switch = campspp_circuit_dynamics(
          computationalcamps(T, N, ones), χ_campspp, thl_campspp, Nmax_campspp, gates, phases, O;
          layer_ends = le, track = true)
        @test s == M
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
        # Frame-dependent: every gate goes to the frame, so bond dimension 1 and k = 0
        # throughout, and CAMPS-PP never switches to PP
        @test s_switch == M && k_switch == 0
        @test all(==(1), bds)
      end
    end
  end
end
