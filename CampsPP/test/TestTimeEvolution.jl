# TimeEvolution.jl against the dense reference at small N, and against exact classical
# results at N = 46. Expected values: CONVENTIONS.md and DenseReference.jl only.
#
# Per-method parameters are kept separate on purpose (χ_camps, χ_campspp, thl_pp, …).

using ITensors, ITensorMPS
import MagicTensors as mt

const bonddim = DisentangleCAMPS.bonddim

# STATE_TYPES, `stackname`, `state_expectation` and `frame_dagger`: helpers.jl.

# -- Test circuits and observables ------------------------------------------------------

"Doped π/4 XXZ circuit as in the production runs: `:xy` doping at 3π/16 on a seeded mask."
function doped_xxz(N, t; p = 0.5, seed = 11, magicphase = 3π/16, kind = :xy)
  gates, phases = xxz_circuit(π/4, π/4, t, N)
  mask = doping_mask(Xoshiro(seed), phases, p)
  phases = kind == :xy ? xy_magic(phases, mask; magicphase) : z_magic(phases, mask; magicphase)
  return gates, phases, mask
end

basis_mps(N, ones) = MPS(ComplexF64, siteinds("Qubit", N), [j in ones ? "1" : "0" for j in 1:N])

"""
Transferred magnetization plus the spin current across the cut,
½(X_a Y_{a+1} − Y_a X_{a+1}) with a = ⌊N/2⌋.

XX, YY, ZZ and the transferred magnetization are real matrices, so on a basis state the
magnetization alone cannot tell ϕ from −ϕ (the state is just complex-conjugated). The
current contains one Y per term and changes sign, which pins the sign of every angle.
"""
function test_observable(N, ones)
  a = fld(N, 2)
  O = transferredmagnetization(N, ones)
  O += pp.PauliString(N, [:X, :Y], [a, a + 1], 0.5)
  O += pp.PauliString(N, [:Y, :X], [a, a + 1], -0.5)
  return O
end

"""
Exact ⟨O⟩ per layer for the undoped π/4 XXZ circuit: every XX-YY-ZZ triple is a SWAP
(CONVENTIONS.md §2; checked on the dense reference in TestCircuits.jl), so the state stays
a basis state whose bits are permuted. ⟨O⟩ = (ones gained by L) + (zeros gained by R).
"""
function swap_dynamics(N, ones, t)
  bits = [j in ones for j in 1:N]
  L = 1:fld(N, 2); R = fld(N, 2)+1:N
  n1L0 = count(bits[L]); n0R0 = count(.!bits[R])
  evs = [0.0]
  for _ in 1:t
    for (a, b) in spec_bonds(N)
      bits[a], bits[b] = bits[b], bits[a]
    end
    push!(evs, (count(bits[L]) - n1L0) + (count(.!bits[R]) - n0R0))
  end
  return evs
end

# The small configurations every method is checked on.
function small_configs()
  cfgs = []
  # A: doped π/4 XXZ (Clifford π/4 gates plus general 3π/16 gates), layer by layer.
  # N = 5 keeps the untruncated Pauli sums of the PP stage small (at most 4^5 terms).
  N, t, ones = 5, 4, [2, 4, 5]
  gates, phases, mask = doped_xxz(N, t)
  push!(cfgs, (name = "doped π/4 XXZ, N=5", N, ones, gates, phases, layer_ends = spec_layer_ends(N, t)))
  # B: XXZ with general angles everywhere, odd N.
  N, t, ones = 5, 3, [1, 4]
  gates, phases = xxz_circuit(0.3, 0.7, t, N)
  push!(cfgs, (name = "general XXZ, N=5", N, ones, gates, phases, layer_ends = spec_layer_ends(N, t)))
  # B′: the same circuit recorded at caller-supplied gate indices that fall inside layers.
  push!(cfgs, (name = "general XXZ, N=5, mid-layer records", N, ones, gates, phases,
               layer_ends = [4, 11, 20, length(gates)]))
  # C: random high-weight rotations, recorded after every gate (no layer ends).
  N, ones = 4, [3]
  gates, phases = rotation_circuit(Xoshiro(21), 10, N)
  @assert all(abs(ϕ - c) > 1e-3 for ϕ in phases, c in (0:8) .* (π / 4))
  push!(cfgs, (name = "random rotations, N=4", N, ones, gates, phases, layer_ends = nothing))
  # D: N = 4, where the last bond of every layer, (2,3), crosses the cut, so a value
  # recorded a few gates early or late is visible.
  N, t, ones = 4, 3, [3]
  gates, phases = xxz_circuit(0.45, -0.8, t, N)
  push!(cfgs, (name = "general XXZ, N=4", N, ones, gates, phases, layer_ends = spec_layer_ends(N, t)))
  return cfgs
end

record_points(cfg) = isnothing(cfg.layer_ends) ? collect(eachindex(cfg.gates)) : cfg.layer_ends

# -- campssrc's layer-end disentangler, rebuilt from the public API ---------------------

"""
The layer-end disentangler of `campssrc_circuit_dynamics` as the old code calls it
(stabSandwich-old, CampsPP/src/TimeEvolution.jl: `cmps.disentangle!(ψ, strat, N;
criterion = crit, min_diff = 6)`, no SVD keywords): up to `max_iter` (there N) sweeps of
`strategy` with `criterion`, stopping after the first sweep that gains ≤ 6. On QubitCAMPS,
the documented equivalent `sweep_disentangler(criterion, strategy, max_iter; min_diff = 6)`.
Returns the number of sweeps run.
"""
function layer_end_disentangle!(ψ::mt.QubitCAMPS, crit, strat, max_iter)
  _, _, sweeps = mt.disentangle!(sweep_disentangler(crit, strat, max_iter; min_diff = 6), ψ)
  return length(sweeps)
end

"""
`campssrc_circuit_dynamics` rebuilt gate by gate on the same stack, for a run where nothing
is truncated: `apply!(ψ, P, ϕ, thl)` per gate and `layer_end_disentangle!` at every record
point in `rec`. Per record point: the frame, the bond dimension, whether the disentangler
changed the frame, the number of sweeps. Also returns the largest |norm² − 1| after any gate.
"""
function src_reconstruction(T, N, ones, gates, phases, thl, rec, crit, strat; max_iter = N)
  ψ = computationalcamps(T, N, ones)
  out = []
  drift = 0.0
  for (i, (g, ϕ)) in enumerate(zip(gates, phases))
    apply!(ψ, PauliOperator(getpauli(g, N)), ϕ, thl)
    drift = max(drift, abs(norm(ψ)^2 - 1))
    if i in rec
      before = frame_dagger(ψ)
      nsweeps = layer_end_disentangle!(ψ, crit, strat, max_iter)
      push!(out, (frame = frame_dagger(ψ), bd = bonddim(ψ), lds = linkdims(mt.get_mps(ψ)),
                  changed = frame_dagger(ψ) != before, nsweeps))
    end
  end
  return out, drift
end

"""
campssrc rebuilt gate by gate with the public API (as `src_reconstruction`), for the χmax
tests: per record point, the bond dimension just before the layer-end disentangler, after
it, and the peak over the gates since the previous record point (that record point's gate
included, before its disentangler).
"""
function src_trace(T, N, ones, gates, phases, thl, rec, crit, strat)
  ψ = computationalcamps(T, N, ones)
  out = NamedTuple{(:before, :after, :peak), NTuple{3, Int}}[]
  peak = 0
  for (i, (g, ϕ)) in enumerate(zip(gates, phases))
    apply!(ψ, PauliOperator(getpauli(g, N)), ϕ, thl)
    peak = max(peak, bonddim(ψ))
    if i in rec
      before = bonddim(ψ)
      layer_end_disentangle!(ψ, crit, strat, N)
      push!(out, (before = before, after = bonddim(ψ), peak = peak))
      peak = 0
    end
  end
  return out
end

@testset "TimeEvolution" begin
  @testset "the doped configuration is really affected by the doping" begin
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = transferredmagnetization(N, ones)
    doped = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    @test maxabsdiff(doped, swap_dynamics(N, ones, length(layer_ends))) > 0.1
    undoped = fill(π/4, length(phases))
    @test maxabsdiff(dense_dynamics(N, ones, gates, undoped, O, layer_ends),
                     swap_dynamics(N, ones, length(layer_ends))) <= EXACT_ATOL
  end

  # CAMPS-PP hand-overs that happen after the first record point, per state type
  late_switch = Dict(T => Bool[] for T in STATE_TYPES)
  early_visible = Bool[] # would a value recorded two gates early differ from the reference?
  late_visible = Bool[]  # … and one gate late?
  # Did the layer-end disentangler change the frame? Per (state type, criterion, strategy).
  SRC_SETTINGS = [(c, s) for c in (:chi3, :entangle) for s in (:snake, :brickwork)]
  frame_changed = Dict((T, c, s) => Bool[] for T in STATE_TYPES for (c, s) in SRC_SETTINGS)

  @testset "exact regime vs dense: $(cfg.name)" for cfg in small_configs()
    (; N, ones, gates, phases, layer_ends) = cfg
    M = length(gates)
    O = test_observable(N, ones)
    rec = record_points(cfg)
    ref = dense_dynamics(N, ones, gates, phases, O, rec)
    @test length(ref) == length(rec) + 1
    @test ref[1] ≈ 0 atol = EXACT_ATOL
    @test maximum(abs, ref) > 0.1                 # the dynamics is not trivial
    # The observable pins the sign of the angles, the gate order and the record points:
    @test maxabsdiff(ref, dense_dynamics(N, ones, gates, -phases, O, rec)) > 0.05
    @test maxabsdiff(ref, dense_dynamics(N, ones, reverse(gates), reverse(phases), O, rec)) > 0.05
    if !isnothing(layer_ends)
      # Sensitivity to a misplaced record point. (One gate early is never visible when a
      # layer ends with ZZ, which commutes with every term of O; two gates early can be.)
      push!(early_visible, maxabsdiff(ref, dense_dynamics(N, ones, gates, phases, O, rec .- 2)) > 1e-3)
      push!(late_visible, maxabsdiff(ref, dense_dynamics(N, ones, gates, phases, O, min.(rec .+ 1, M))) > 1e-3)
    end
    @testset "camps_circuit_dynamics, $(stackname(T))" for T in STATE_TYPES
      ψ0 = computationalcamps(T, N, ones)
      χ_camps = 10^6                              # never reached: runs the whole circuit
      ψ, k, s, evs = camps_circuit_dynamics(ψ0, gates, phases, χ_camps, O; layer_ends)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test 0 <= k <= N
      ev = state_expectation(ψ, O)                # CONVENTIONS.md §5: a real number
      @test abs(imag(ev)) <= EXACT_ATOL
      @test real(ev) ≈ ref[end] atol = EXACT_ATOL
      # The `track` form returns the same values, plus one bond dimension per record.
      ψt, kt, st, evst, bds = camps_circuit_dynamics(ψ0, gates, phases, χ_camps, O; layer_ends, track = true)
      @test evst == evs && kt == k && st == s
      @test length(bds) == length(evs) && bds[1] == 1 && bds[end] == bonddim(ψt)
    end

    @testset "pauliprop_circuit_dynamics (PP alone)" begin
      thl_pp, Nmax_pp = 0.0, typemax(Int)         # no truncation
      s, evs = pauliprop_circuit_dynamics(ones, gates, phases, thl_pp, Nmax_pp, O; layer_ends)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
    end

    @testset "campspp_circuit_dynamics, χ = $χ_campspp, $(stackname(T))" for T in STATE_TYPES,
        χ_campspp in (2, 4, 10^6)
      ψ0 = computationalcamps(T, N, ones)
      thl_campspp, Nmax_campspp = 0.0, typemax(Int)
      evs, s, bonddims, nterms, _, _, s_switch, k_switch = campspp_circuit_dynamics(
        ψ0, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, O; layer_ends, track = true)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL    # exact on both sides of the switch
      # The non-`track` form returns the same values.
      evs2, s2 = campspp_circuit_dynamics(
        ψ0, χ_campspp, thl_campspp, Nmax_campspp, gates, phases, O; layer_ends)
      @test evs2 == evs && s2 == s

      # Frame-dependent: where the hand-over happens. It is the first gate after which the
      # CAMPS bond dimension is ≥ χ: reached after s_switch gates, not after s_switch − 1.
      ψc, kc, sc, _ = camps_circuit_dynamics(ψ0, gates, phases, χ_campspp, O; layer_ends)
      @test sc == s_switch && kc == k_switch
      if χ_campspp == 10^6
        @test s_switch == M                       # never switched
      else
        χ_campspp == 2 && @test 0 < s_switch < M  # every configuration switches at χ = 2
        if s_switch < M
          @test bonddim(ψc) >= χ_campspp
          @test bonddim(evolve(computationalcamps(T, N, ones), s_switch - 1, gates, phases)[1]) < χ_campspp
          push!(late_switch[T], s_switch > first(rec))
        end
      end
    end

    @testset "campssrc_circuit_dynamics ($crit, $strat), $(stackname(T))" for T in STATE_TYPES,
        (crit, strat) in SRC_SETTINGS
      # SVD cutoff. Nothing is discarded here: below, the norm² moves by rounding only
      # (≤ 1e-12, far below thl_camps). Not 1e-15: campssrc stops disentangling once the
      # norm² has dropped by more than a fraction thl, and at thl ≲ 1e-13 rounding alone
      # (~1e-14) can trip that before the first layer end, so ⟨O⟩, which is
      # frame-invariant, would pass with a broken or skipped layer-end disentangler (T9
      # review). What rounding does at such thresholds is not specified, and not tested.
      thl_camps = 1e-11
      ψ0T = computationalcamps(T, N, ones)
      ψ, s, evs, bds = campssrc_circuit_dynamics(ψ0T, gates, phases, thl_camps, O;
        criterion = crit, strategy = strat, layer_ends, track = true)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test abs(imag(state_expectation(ψ, O))) <= EXACT_ATOL
      @test length(bds) == length(evs) && bds[1] == 1 && bds[end] == bonddim(ψ)
      # Frame-dependent: the bond dimension at every record point and the final frame equal
      # those of the same circuit rebuilt with the public API, the old code's disentangler
      # run at every record point (`src_reconstruction`, same stack: with :entangle, ties are
      # broken by rounding, so frames are compared within a stack only). This catches a
      # disentangler that is skipped, or set up with another strategy, criterion or stop
      # rule.
      recon, drift = src_reconstruction(T, N, ones, gates, phases, thl_camps, rec, crit, strat)
      @test drift <= 1e-12                        # rounding only: nothing was truncated
      @test bds == [1; [r.bd for r in recon]]
      @test frame_dagger(ψ) == recon[end].frame
      if crit == :entangle
        # The frame at every record point: that of the run cut there. (Only where the
        # disentangler acts: each cut run repeats the disentangler from the start, so this
        # costs O(record points²).)
        src_frames = map(eachindex(rec)) do l
          le_l = isnothing(layer_ends) ? nothing : layer_ends[1:l]
          ψl, = campssrc_circuit_dynamics(ψ0T, gates[1:rec[l]], phases[1:rec[l]], thl_camps, O;
            criterion = crit, strategy = strat, layer_ends = le_l)
          frame_dagger(ψl)
        end
        @test src_frames == [r.frame for r in recon]
      end
      # Not vacuous: with :entangle the disentangler changes the frame at (almost) every
      # record point, so skipping it anywhere would show. (On QubitCAMPS, one layer end of
      # D with :brickwork finds nothing to gain; allow one such point per configuration.)
      crit == :entangle && @test count(r.changed for r in recon) >= length(rec) - 1
      # Frame-dependent: compare the frame with the same gates applied without the
      # layer-end disentangler, and the bond dimension at every record point.
      plain = computationalcamps(T, N, ones)
      plain_bds = [bonddim(plain)]
      plain_lds = Vector{Int}[]                   # every bond, at each record point
      for (i, (g, ϕ)) in enumerate(zip(gates, phases))
        apply!(plain, PauliOperator(getpauli(g, N)), ϕ, thl_camps)
        if i in rec
          push!(plain_bds, bonddim(plain))
          push!(plain_lds, linkdims(mt.get_mps(plain)))
        end
      end
      push!(frame_changed[(T, crit, strat)], frame_dagger(ψ) != frame_dagger(plain))
      # [item 8] In the doped configuration, the layer-end bond dimensions with :chi3 are at
      # most those of the run without the disentangler (plan item 8, "Dynamics"), bond by
      # bond at every record point (the per-bond dimensions are those of `recon`, which
      # matches campssrc above), and strictly lower on some bond at some record point, so
      # the check sees a reduction. (On this circuit: bond 2 at the first layer end, 2 vs 4.
      # The maximum alone, `bds`, is non-strict and could not tell.)
      if crit == :chi3 && startswith(cfg.name, "doped")
        @test all(bds .<= plain_bds)
        @test length(plain_lds) == length(recon)
        @test all(all(r.lds .<= p) for (r, p) in zip(recon, plain_lds))
        @test any(any(r.lds .< p) for (r, p) in zip(recon, plain_lds))
      end
    end

    @testset "mps_circuit_dynamics" begin
      thl_mps = 1e-15
      ψ, evs = mps_circuit_dynamics(basis_mps(N, ones), gates, phases, thl_mps, O, nothing; layer_ends)
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
    end
  end

  @testset "the small configurations exercise the hand-over and the disentangler" begin
    # Guards against vacuous passes above (frame-dependent facts, checked per stack).
    @test count(early_visible) >= 2     # record points two gates early would be noticed
    @test count(late_visible) >= 2      # … and one gate late
    for T in STATE_TYPES
      @test count(late_switch[T]) >= 2  # hand-overs after the first record point
      # The layer-end disentangler does change the frame with :entangle.
      @test count(frame_changed[(T, :entangle, :brickwork)]) >= 2
      @test count(frame_changed[(T, :entangle, :snake)]) >= 2
      # [item 8] With :chi3 ((:chi3, :snake) is the production setting; min_diff = 6) it now
      # acts too: singular values ≤ τ·σ_max count as zero (CONVENTIONS.md §2). In the doped
      # configuration (the first of small_configs()) it changes the frame with both
      # strategies. This replaces the old no-op guard `!any(frame_changed[(T, :chi3, …)])`.
      # Tested directly in DisentangleCAMPS/TestDisentangleQubitCAMPS.jl.
      @test startswith(first(small_configs()).name, "doped")
      @test frame_changed[(T, :chi3, :snake)][1]
      @test frame_changed[(T, :chi3, :brickwork)][1]
    end
  end

  @testset "campssrc: the layer-end disentangler runs up to N sweeps, $(stackname(T))" for T in STATE_TYPES
    # The old code calls the disentangler with max_iter = N and min_diff = 6, so a second
    # sweep runs only after a sweep that gains more than 6 (for :entangle, bits of bond
    # entropy summed over the bonds). At N ≤ 5 that cannot happen: the bond entropies sum
    # to at most Σ_b min(b, N − b) = 1 + 2 + 2 + 1 = 6. So the small configurations cannot
    # tell max_iter = N from max_iter = 1; this one can.
    # N = 12, four brickwork layers of exp(i(π/4 + ε)·P_b Q_{b+1}) with random P, Q ∈
    # {X, Y, Z}: each gate is a Clifford entangler times exp(iε·PQ), so (ϕ is not a Clifford
    # angle) the MPS part holds much entanglement that Cliffords can remove. One record
    # point, at the end. ε = 0.05 is ≥ 1e-3 away from every multiple of π/4 (§2).
    N, depth, ε = 12, 4, 0.05
    rng = Xoshiro(1)
    gates, phases = pp.PauliRotation[], Float64[]
    for d in 1:depth, b in (isodd(d) ? (1:2:N-1) : (2:2:N-1))
      push!(gates, pp.PauliRotation([rand(rng, [:X, :Y, :Z]), rand(rng, [:X, :Y, :Z])], [b, b + 1]))
      push!(phases, π/4 + ε)
    end
    ones = findall(rand(rng, Bool, N))
    M = length(gates)
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, [M])
    thl_camps = 1e-11                             # nothing truncated: see the drift check
    for strat in (:brickwork, :snake)
      ψ, s, evs, bds = campssrc_circuit_dynamics(computationalcamps(T, N, ones), gates, phases,
        thl_camps, O; criterion = :entangle, strategy = strat, layer_ends = [M], track = true)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      recon, drift = src_reconstruction(T, N, ones, gates, phases, thl_camps, [M], :entangle, strat)
      @test drift <= 1e-12                        # rounding only
      # Not vacuous: the first sweep gains > 6, a second one runs and changes the frame.
      @test recon[1].nsweeps >= 2
      one_sweep, _ = src_reconstruction(T, N, ones, gates, phases, thl_camps, [M], :entangle,
        strat; max_iter = 1)
      @test one_sweep[1].frame != recon[1].frame
      # Frame-dependent: campssrc ran the sweeps the old code runs (same stack).
      @test frame_dagger(ψ) == recon[1].frame
      @test bds == [1, recon[1].bd]
    end
  end

  @testset "CAMPS-PP continues from a CAMPS state, layer by layer, $(stackname(T))" for T in STATE_TYPES
    # Evolve s0 gates with CAMPS, hand over to PP for the rest. PP must record the dense
    # value at every later layer end, wherever the hand-over is. (No t = 0 entry here.)
    # N = 4 keeps the untruncated Pauli sums small (at most 4^4 terms).
    N, t, ones = 4, 3, [3, 4]
    gates, phases, _ = doped_xxz(N, t)
    M = length(gates)
    le = spec_layer_ends(N, t)
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, le)          # t = 0, then the 3 layers
    @test maxabsdiff(ref, dense_dynamics(N, ones, gates, -phases, O, le)) > 0.05
    for s0 in (0, 5, le[1], le[1] + 1, 14, le[2], M - 1, M)
      ψ, k = evolve(computationalcamps(T, N, ones), s0, gates, phases)
      s, evs = pauliprop_circuit_dynamics(ψ, s0, gates, phases, 0.0, typemax(Int), O; layer_ends = le)
      later = [l for l in 1:t if le[l] > s0]                      # layers still to come
      @test s == M
      @test length(evs) == length(later)
      isempty(later) || @test maxabsdiff(evs, ref[later .+ 1]) <= EXACT_ATOL
    end
  end

  @testset "CAMPS-PP: `thl` and `Nmax` reach the PP stage, $(stackname(T))" for T in STATE_TYPES
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    M = length(gates)
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    ψ0 = computationalcamps(T, N, ones)
    χ_campspp = 2

    # Small Nmax: PP stops at the first layer end where the sum has ≥ Nmax terms. What was
    # recorded before is exact.
    # (The observable itself has 8 terms, so Nmax must be larger than that to be informative.)
    Nmax_campspp = 50
    evs, s, _, nterms, _, _, s_switch, _ = campspp_circuit_dynamics(
      ψ0, χ_campspp, 0.0, Nmax_campspp, gates, phases, O; layer_ends, track = true)
    @test s_switch < layer_ends[1]                # PP sees both layer ends below
    @test s == layer_ends[2]                      # below Nmax at layer 1, above at layer 2
    @test length(evs) == 3
    @test length(nterms) == 3 && nterms[1] == length(O)
    @test nterms[2] < Nmax_campspp <= nterms[3]
    @test maxabsdiff(evs, ref[1:length(evs)]) <= EXACT_ATOL

    # Loose threshold: runs to the end, with a visible truncation error after the switch.
    evs, s, _, _, _, _, s_switch, _ = campspp_circuit_dynamics(
      ψ0, χ_campspp, 0.3, typemax(Int), gates, phases, O; layer_ends, track = true)
    @test s == M && s_switch < M
    @test length(evs) == length(ref)
    @test maxabsdiff(evs, ref) > 1e-3
  end

  @testset "truncated regime: PP threshold" begin
    # A truncated run must approach the dense result as the threshold is tightened.
    # No absolute bound is asserted for the loose thresholds (none is specified, §6).
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    runs = map((3e-1, 3e-2, 1e-3, 1e-6, 0.0)) do thl_pp
      pauliprop_circuit_dynamics(ones, gates, phases, thl_pp, typemax(Int), O; layer_ends)[2]
    end
    @test all(length(evs) == length(ref) for evs in runs)
    errs = [maxabsdiff(evs, ref) for evs in runs]
    @test errs[end] <= EXACT_ATOL
    @test errs[1] > 1e-3                          # the loosest threshold really truncates
    @test nonincreasing(errs)              # error does not grow as thl shrinks
  end

  @testset "truncated regime: CAMPS SVD threshold, $(stackname(T))" for T in STATE_TYPES
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    ψ0 = computationalcamps(T, N, ones)
    runs = map((1e-1, 3e-2, 1e-2, 3e-3, 1e-15)) do thl_camps
      campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
        criterion = :chi3, strategy = :snake, layer_ends)[3]
    end
    @test all(length(evs) == length(ref) for evs in runs)
    errs = [maxabsdiff(evs, ref) for evs in runs]
    @test errs[end] <= EXACT_ATOL
    @test errs[1] > 1e-3
    @test nonincreasing(errs)
    # More than one threshold truncates, by different amounts. (The expectation is not
    # normalised, so part of the error is the norm lost to the cutoff.)
    @test length(unique(round.(errs[errs .> 1e-6]; sigdigits = 3))) >= 2
  end

  @testset "truncated regime: MPS SVD threshold" begin
    # mps_circuit_dynamics truncates after every gate at `cutoff = thl` and renormalises, as
    # the old code (`ITensorMPS.apply(ψ, gate; cutoff = thl); normalize!(ψ)`), so its ⟨O⟩
    # is that of a normalised state. Configuration A (N = 5).
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    runs = map((3e-1, 1e-1, 3e-2, 1e-3, 1e-15)) do thl_mps
      mps_circuit_dynamics(basis_mps(N, ones), gates, phases, thl_mps, O, nothing; layer_ends)
    end
    @test all(length(evs) == length(ref) for (_, evs) in runs)
    errs = [maxabsdiff(evs, ref) for (_, evs) in runs]
    @test errs[end] <= EXACT_ATOL                 # vanishes without truncation
    @test errs[1] > 1e-3                          # the loosest threshold really truncates
    @test nonincreasing(errs)                     # error does not grow as thl shrinks
    # At least two thresholds really truncate: the final bond dimension is below that of
    # the untruncated run (whose ⟨O⟩ is exact, above).
    χs = [maxlinkdim(ψ) for (ψ, _) in runs]
    @test count(<(χs[end]), χs) >= 2
    # The returned state is normalised, also after truncation (rounding only).
    @test all(abs(norm(ψ) - 1) <= 1e-12 for (ψ, _) in runs)
  end

  @testset "truncated regime: Nmax stops PP early" begin
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    le = layer_ends
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, le)
    # Term counts of the untruncated run choose an Nmax that is first reached at layer 2.
    _, _, full = pauliprop_circuit_dynamics(ones, gates, phases, 0.0, typemax(Int), O;
      layer_ends = le, track = true)
    @test length(full) == length(le) + 1
    @test full[2] < full[3]
    Nmax_pp = full[3]
    s, evs, nterms = pauliprop_circuit_dynamics(ones, gates, phases, 0.0, Nmax_pp, O;
      layer_ends = le, track = true)
    @test s == le[2]                              # stopped at the end of layer 2
    @test length(evs) == 3 == length(nterms)      # t = 0, layer 1, layer 2
    @test nterms == full[1:3]
    @test maxabsdiff(evs, ref[1:3]) <= EXACT_ATOL # what was computed is exact
  end

  @testset "output file, $(stackname(T))" for T in STATE_TYPES
    # With a file name as `output`, every recorded value is appended as "index<TAB>value".
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    ψ0 = computationalcamps(T, N, ones)
    mktempdir() do dir
      path = joinpath(dir, "evs.txt")
      evs, s = campspp_circuit_dynamics(ψ0, 2, 0.0, typemax(Int), gates, phases, O, path; layer_ends)
      rows = [split(l, '\t') for l in eachline(path) if !isempty(strip(l)) && !startswith(l, '#')]
      @test [parse(Int, r[1]) for r in rows] == collect(0:length(layer_ends))
      @test maxabsdiff([parse(Float64, r[2]) for r in rows], evs) <= 1e-12
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test count(startswith("# CAMPS stopped"), eachline(path)) == 1
      @test count(startswith("# Pauli prop. stopped"), eachline(path)) == 1
    end
  end

  @testset "`:z` doping acts trivially at ϕ_xy = ϕ_z = π/4" begin
    # A doped triple is SWAP·exp(i(α−π/4)ZZ): on a basis state the ZZ factor is a phase, so
    # every ⟨O⟩ equals the undoped value, for any mask and any angle α.
    @testset "dense, N=6" begin
      N, t, ones = 6, 4, [2, 4, 5, 6]
      O = transferredmagnetization(N, ones)
      le = spec_layer_ends(N, t)
      exact = swap_dynamics(N, ones, t)
      for α in (π/3, 3π/16, 1.234)
        gates, phases, mask = doped_xxz(N, t; p = 0.5, kind = :z, magicphase = α)
        @assert any(mask)
        @test maxabsdiff(dense_dynamics(N, ones, gates, phases, O, le), exact) <= EXACT_ATOL
      end
    end

    @testset "code, N=46, $t layers" for t in (get(ENV, "STABSANDWICH_SLOW", "") == "" ? (3,) : (3, 23))
      N = 46
      _, ones = domainwallstate(Xoshiro(100), N, 0.6)
      O = transferredmagnetization(N, ones)
      le = spec_layer_ends(N, t)
      exact = swap_dynamics(N, ones, t)
      gates, undoped = xxz_circuit(π/4, π/4, t, N)
      M = length(gates)
      χ_campspp, thl_campspp, Nmax_campspp = 64, 1e-10, 10^6
      thl_camps = 1e-12

      run_pp(T, phases) = campspp_circuit_dynamics(computationalcamps(T, N, ones), χ_campspp,
        thl_campspp, Nmax_campspp, gates, phases, O; layer_ends = le, track = true)
      run_src(T, phases) = campssrc_circuit_dynamics(computationalcamps(T, N, ones), gates, phases,
        thl_camps, O; criterion = :chi3, strategy = :snake, layer_ends = le, track = true)

      pp0 = Dict(T => run_pp(T, undoped)[1] for T in STATE_TYPES)
      # undoped runs against the exact SWAP result
      @test all(maxabsdiff(pp0[T], exact) <= EXACT_ATOL for T in STATE_TYPES)
      src0 = Dict(T => run_src(T, undoped)[3] for T in STATE_TYPES)
      @test all(maxabsdiff(src0[T], exact) <= EXACT_ATOL for T in STATE_TYPES)

      for α in (π/3, 1.234)
        _, phases, mask = doped_xxz(N, t; p = 0.3, kind = :z, magicphase = α)
        @assert count(mask) > 10
        @testset "campspp, $(stackname(T))" for T in STATE_TYPES
          evs, s, bds, _, _, _, s_switch, k = run_pp(T, phases)
          @test evs == pp0[T]                      # bitwise equal to the p = 0 run
          @test s == M && s_switch == M            # no switch to PP
          @test all(==(1), bds) && k == 0          # frame-dependent: nothing was added
        end
        # The truncating method sends every doped gate through an SVD and never
        # renormalises, and `expectation` is not normalised (CONVENTIONS.md §5). Over a few
        # hundred doped gates the norm drifts by rounding (measured 4e-12 on norm² at 23
        # layers), which scales ⟨O⟩ ≈ 30 by the same factor. So: the norm stays within 1e-10
        # of 1, each layer agrees up to that drift, and the normalised final value agrees.
        @testset "campssrc, $(stackname(T))" for T in STATE_TYPES
          ψsrc, _, src, srcbds = run_src(T, phases)
          drift = abs(norm(ψsrc)^2 - 1)
          @test drift <= EXACT_ATOL
          @test all(abs.(src .- src0[T]) .<= EXACT_ATOL .+ abs.(src0[T]) .* drift)
          @test abs(src[end] / norm(ψsrc)^2 - src0[T][end]) <= EXACT_ATOL
          @test all(==(1), srcbds)                 # frame-dependent
        end
      end
    end
  end

  @testset "large N, Clifford only: undoped π/4 XXZ at N=46, $(stackname(T))" for T in STATE_TYPES
    N, t = 46, get(ENV, "STABSANDWICH_SLOW", "") == "" ? 6 : 23
    _, ones = domainwallstate(Xoshiro(100), N, 0.6)
    O = transferredmagnetization(N, ones)
    le = spec_layer_ends(N, t)
    exact = swap_dynamics(N, ones, t)
    @test maximum(abs, exact) >= 2                # excitations do cross the cut
    gates, phases = xxz_circuit(π/4, π/4, t, N)
    ψ0 = computationalcamps(T, N, ones)

    ψ, k, s, evs = camps_circuit_dynamics(ψ0, gates, phases, 64, O; layer_ends = le)
    @test s == length(gates) && k == 0
    @test maxabsdiff(evs, exact) <= EXACT_ATOL
    @test DisentangleCAMPS.bonddim(ψ) == 1         # frame-dependent

    if T === first(STATE_TYPES)                    # PP alone has no CAMPS state: run once
      s, evs = pauliprop_circuit_dynamics(ones, gates, phases, 1e-10, 10^6, O; layer_ends = le)
      @test s == length(gates)
      @test maxabsdiff(evs, exact) <= EXACT_ATOL
    end
  end

  # -- Kicked Ising: every driver against the dense reference (KIM plan step 2) -----------
  # N = 6, t = 6, from the Néel state, observable local_z(N, j) for every j; untruncated
  # (thresholds 0, or ≤ 1e-12 where a threshold drives other logic), so exact to 1e-10
  # (§6). The reference is the spec circuit (KimSpec.jl) on the dense simulator; the
  # drivers get `kicked_ising_circuit`. ε: 0 (Clifford), a seeded uniform value, and a
  # seeded per-site vector with one Clifford site (ε_3 = 0). Non-zero ε are ≥ 2e-3 in size,
  # so the Z angles are ≥ 1e-3 from every multiple of π/4 (§2, tests' angles).
  # t = 6, not 4: ⟨Z_j⟩ is 0 after layers 1 … N/2 = 3 for any ε (TestCircuits.jl), so the
  # comparison sees the circuit at layers 4–6 only.
  # At ε = 0 every recorded ⟨Z_j⟩ after t = 0 is 0 (the zero window lasts to layer 11 at
  # N = 6), and so is ⟨Z_j⟩ from many wrong Clifford circuits: that configuration checks the
  # frame (no switch, k = 0), not values. The Clifford point's values are checked at N = 46
  # with a sensitive observable (TestCliffordLargeN.jl).
  kim_N, kim_t = 6, 6
  kim_εs = let εv = 0.5 .* (2 .* rand(Xoshiro(32), kim_N) .- 1)
    εv[3] = 0.0
    [("ε = 0 (Clifford)", 0.0), ("uniform ε", 0.1 + 0.4 * rand(Xoshiro(31))), ("per-site ε", εv)]
  end
  kim_profiles = Dict(name => kim_dense_profile(ε, kim_t, kim_N) for (name, ε) in kim_εs)   # (t + 1) × N

  @testset "KIM driver references, N=6: $name" for (name, ε) in kim_εs
    N, t = kim_N, kim_t
    prof = kim_profiles[name]
    @test all(e -> e == 0 || abs(e) >= 2e-3, ε)
    @test prof[1, :] == [j in kim_neel(N) ? -1.0 : 1.0 for j in 1:N]
    if ε != 0
      # The reference is not 0 after t = 0, sees ε, and ε on the wrong qubit (per-site case)
      @test maximum(abs, prof[2:end, :]) > 0.05
      @test maximum(abs.(prof .- kim_dense_profile(0.0, t, N))) > 0.05
      ε isa AbstractVector &&
        @test maximum(abs.(prof .- kim_dense_profile(circshift(ε, 1), t, N))) > 1e-3
    end
  end

  kim_switched = Bool[]                          # CAMPS-PP at χ = 4 switched mid-circuit

  @testset "KIM drivers vs dense, N=6: $name" for (name, ε) in kim_εs
    N, t = kim_N, kim_t
    ones = kim_neel(N)
    gates, phases = kicked_ising_circuit(ε, t, N)
    M = length(gates)
    layer_ends = [kim_spec_L(N) * l for l in 1:t]
    prof = kim_profiles[name]

    @testset "j = $j" for j in 1:N
      O = local_z(N, j)
      ref = prof[:, j]

      @testset "campspp_circuit_dynamics, χ = $χ_campspp, $(stackname(T))" for T in STATE_TYPES,
          χ_campspp in (4, 10^6)
        thl_campspp, Nmax_campspp = 0.0, typemax(Int)
        evs, s, _, _, _, _, s_switch, k_switch = campspp_circuit_dynamics(
          computationalcamps(T, N, ones), χ_campspp, thl_campspp, Nmax_campspp, gates, phases, O;
          layer_ends, track = true)
        @test s == M
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
        # Frame-dependent: at ε = 0 every gate is Clifford, so the bond dimension stays 1
        # and k = 0: no switch at any χ > 1
        ε == 0 && @test s_switch == M && k_switch == 0
        χ_campspp == 10^6 && @test s_switch == M
        χ_campspp == 4 && push!(kim_switched, 0 < s_switch < M)
      end

      @testset "camps_circuit_dynamics, $(stackname(T))" for T in STATE_TYPES
        χ_camps = 10^6
        _, _, s, evs = camps_circuit_dynamics(computationalcamps(T, N, ones), gates, phases,
                                               χ_camps, O; layer_ends)
        @test s == M
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
      end

      @testset "campssrc_circuit_dynamics ($crit, $strat), $(stackname(T))" for T in STATE_TYPES,
          (crit, strat) in [(:chi3, :snake), (:entangle, :brickwork)]
        thl_camps = 1e-12
        _, s, evs = campssrc_circuit_dynamics(computationalcamps(T, N, ones), gates, phases,
          thl_camps, O; criterion = crit, strategy = strat, layer_ends)
        @test s == M
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
      end

      @testset "mps_circuit_dynamics" begin
        thl_mps = 1e-15
        _, evs = mps_circuit_dynamics(basis_mps(N, ones), gates, phases, thl_mps, O, nothing;
                                      layer_ends)
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
      end

      @testset "pauliprop_circuit_dynamics" begin
        thl_pp, Nmax_pp = 0.0, typemax(Int)
        s, evs = pauliprop_circuit_dynamics(ones, gates, phases, thl_pp, Nmax_pp, O; layer_ends)
        @test s == M
        @test maxabsdiff(evs, ref) <= EXACT_ATOL
      end
    end
  end

  @testset "KIM: CAMPS-PP at χ = 4 does switch mid-circuit" begin
    @test count(kim_switched) >= 2              # so both sides of the switch are checked
  end

  # -- Budget revision: the `χmax` keyword (CONVENTIONS.md §2 "Recording") -----------------
  # `campssrc_circuit_dynamics` and `mps_circuit_dynamics` take `χmax`. With `nothing` (the
  # default) a run goes to the end; with an integer it stops after the first record point
  # whose bond dimension (largest link dimension; campssrc: after the layer-end disentangler)
  # is ≥ χmax, and that value is still recorded. A larger bond dimension between record points
  # doesn't stop a run.
  #
  # Expected values: the same driver without χmax (its unchanged, tested behaviour), cut at
  # the record point predicted from bond dimensions taken from other sources: for campssrc,
  # the unstopped run's `track` output, checked against the public-API rebuild `src_trace`
  # (which also gives the bond dimension before the disentangler and the peaks between record
  # points); for MPS, `maxlinkdim` of the states of unstopped runs on the first l layers. χmax
  # takes every recorded value, so a `>` rule, a check before the disentangler or a check
  # after every gate predicts a different stop somewhere (non-vacuity counted below).
  χ_cfgs = Any[]
  # The first two XXZ configs, and one recorded after every gate (`layer_ends = nothing`), so a
  # stop rule that works only with an explicit list of record points fails there
  every_gate = only(filter(c -> c.layer_ends === nothing, small_configs()))
  for cfg in [small_configs()[1:2]; every_gate]
    push!(χ_cfgs, (cfg.name, cfg.N, cfg.ones, cfg.gates, cfg.phases,
                   test_observable(cfg.N, cfg.ones), cfg.layer_ends))
  end
  for (name, ε) in kim_εs[2:3]
    g, ϕ = kicked_ising_circuit(ε, kim_t, kim_N)
    push!(χ_cfgs, ("KIM N=6, $name", kim_N, kim_neel(kim_N), g, ϕ, local_z(kim_N, 1),
                   [kim_spec_L(kim_N) * l for l in 1:kim_t]))
  end
  # A circuit with a bond-dimension peak between record points only: two magic qubits,
  # then XX on them and its inverse. The state after the pair is a product state again.
  peak_gates = [pp.PauliRotation([:X], [1]), pp.PauliRotation([:X], [2]),
                pp.PauliRotation([:X, :X], [1, 2]), pp.PauliRotation([:X, :X], [1, 2])]
  push!(χ_cfgs, ("mid-layer peak, N=3", 3, Int[], repeat(peak_gates, 2),
                 repeat([0.37, 0.29, 0.41, -0.41], 2), transferredmagnetization(3, Int[]),
                 [2, 4, 6, 8]))

  # Non-vacuity: does some χmax tell the spec's stop from a wrong rule's?
  seen = Dict(k => false for k in ("stops early", "≥ vs >", "after vs before disentangler",
                                   "record points vs every gate", "mps stops early",
                                   "mps ≥ vs >", "mps mid-layer peak",
                                   "stops early, layer_ends = nothing",
                                   "mps stops early, layer_ends = nothing"))
  # First record point with value ≥ χ (nothing: none)
  first_ge(vals, χ) = findfirst(>=(χ), vals)

  @testset "χmax, campssrc ($crit, $strat): $name, $(stackname(T))" for T in STATE_TYPES,
      (crit, strat) in [(:chi3, :snake), (:entangle, :brickwork)],
      (name, N, ones, gates, phases, O, rec) in χ_cfgs
    M = length(gates)
    recv = rec === nothing ? collect(1:M) : rec  # §2: nothing = every gate is a record point
    thl_camps = 1e-11                           # ≥ 1e-11: the disentangler runs (§2)
    ψ0 = computationalcamps(T, N, ones)
    ψf, sf, evsf, bdsf = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
      criterion = crit, strategy = strat, layer_ends = rec, track = true)
    @test sf == M
    # χmax = nothing is the call without the keyword
    ψn, sn, evsn, bdsn = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
      criterion = crit, strategy = strat, layer_ends = rec, track = true, χmax = nothing)
    @test sn == sf && evsn == evsf && bdsn == bdsf && frame_dagger(ψn) == frame_dagger(ψf)

    tr = src_trace(T, N, ones, gates, phases, thl_camps, recv, crit, strat)
    after = [r.after for r in tr]
    @test bdsf == [1; after]                    # track output = the rebuild, after the disentangler
    before = [r.before for r in tr]
    peaks = [r.peak for r in tr]
    name == "mid-layer peak, N=3" && @test maximum(peaks) > maximum(after)   # the premise

    for χ in sort(unique([2; after; before; peaks; maximum(after) + 1]))
      χ >= 2 || continue
      l = first_ge(after, χ)
      ψs, ss, evss, bdss = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
        criterion = crit, strategy = strat, layer_ends = rec, track = true, χmax = χ)
      if l === nothing                          # never reached: the whole circuit
        @test ss == M && evss == evsf && bdss == bdsf
      else
        @test ss == recv[l]
        @test length(evss) == l + 1 && maxabsdiff(evss, evsf[1:l+1]) <= 1e-12
        @test bdss == bdsf[1:l+1]
        @test bonddim(ψs) == after[l]
        l < length(recv) && (seen["stops early"] = true)
        l < length(recv) && rec === nothing && (seen["stops early, layer_ends = nothing"] = true)
      end
      first_ge(after .- 1, χ) != l && (seen["≥ vs >"] = true)    # > χ ⟺ ≥ χ + 1
      first_ge(before, χ) != l && (seen["after vs before disentangler"] = true)
      first_ge(peaks, χ) != l && (seen["record points vs every gate"] = true)
    end
  end

  @testset "χmax, mps_circuit_dynamics: $name" for (name, N, ones, gates, phases, O, rec) in χ_cfgs
    M = length(gates)
    recv = rec === nothing ? collect(1:M) : rec
    thl_mps = 1e-15
    mrun(g, ϕ, le; kw...) = mps_circuit_dynamics(basis_mps(N, ones), g, ϕ, thl_mps, O, nothing;
                                                 layer_ends = le, kw...)
    ψf, evsf = mrun(gates, phases, rec)
    ψn, evsn = mrun(gates, phases, rec; χmax = nothing)
    @test evsn == evsf && maxlinkdim(ψn) == maxlinkdim(ψf)
    # Bond dimension at record point l: the state of the unstopped run on gates 1 … recv[l]
    bd = [maxlinkdim(first(mrun(gates[1:recv[l]], phases[1:recv[l]], recv[1:l]))) for l in eachindex(recv)]
    if name == "mid-layer peak, N=3"
      # Premise: after gate 3 (XX) the bond dimension is 2, at every record point 1
      @test maxlinkdim(first(mrun(gates[1:3], phases[1:3], [3]))) == 2
      @test all(==(1), bd)
      ψs, evss = mrun(gates, phases, rec; χmax = 2)
      @test evss == evsf                        # the peak alone doesn't stop
      seen["mps mid-layer peak"] = true
    end
    for χ in sort(unique([2; bd; maximum(bd) + 1]))
      l = first_ge(bd, χ)
      ψs, evss = mrun(gates, phases, rec; χmax = χ)
      if l === nothing
        @test evss == evsf
      else
        @test length(evss) == l + 1 && maxabsdiff(evss, evsf[1:l+1]) <= 1e-12
        @test maxlinkdim(ψs) == bd[l]
        l < length(recv) && (seen["mps stops early"] = true)
        l < length(recv) && rec === nothing && (seen["mps stops early, layer_ends = nothing"] = true)
      end
      first_ge(bd .- 1, χ) != l && (seen["mps ≥ vs >"] = true)
    end
  end

  @testset "χmax: the output file's stop line" begin
    # The stop line names the bond-dimension reason, at the stop gate; without χmax, the end
    # of the circuit (the wording around the reason is the drivers' own, not pinned)
    name, N, ones, gates, phases, O, rec = χ_cfgs[1]
    M = length(gates)
    ψ0 = computationalcamps(N, ones)
    tr = src_trace(mt.QubitCAMPS, N, ones, gates, phases, 1e-11, rec, :chi3, :snake)
    after = [r.after for r in tr]
    # the smallest χmax ≥ 2 that stops before the last record point
    χ = minimum(filter(c -> c >= 2 && something(findfirst(>=(c), after), length(rec)) < length(rec), after))
    l = findfirst(>=(χ), after)
    stoplines(path) = filter(l -> occursin("stopped at gate", l), readlines(path))
    mktempdir() do dir
      p1, p2 = joinpath(dir, "src_stop.txt"), joinpath(dir, "src_end.txt")
      campssrc_circuit_dynamics(ψ0, gates, phases, 1e-11, O, p1; layer_ends = rec, χmax = χ)
      campssrc_circuit_dynamics(ψ0, gates, phases, 1e-11, O, p2; layer_ends = rec)
      @test length(stoplines(p1)) == 1 && occursin("stopped at gate $(rec[l]) ", only(stoplines(p1))) &&
            occursin("bond dim", only(stoplines(p1)))
      @test length(stoplines(p2)) == 1 && occursin("stopped at gate $M ", only(stoplines(p2))) &&
            occursin("end of circuit", only(stoplines(p2)))

      bd = [maxlinkdim(first(mps_circuit_dynamics(basis_mps(N, ones), gates[1:rec[k]], phases[1:rec[k]],
              1e-15, O, nothing; layer_ends = rec[1:k]))) for k in eachindex(rec)]
      χm = minimum(filter(c -> c >= 2 && something(findfirst(>=(c), bd), length(rec)) < length(rec), bd))
      lm = findfirst(>=(χm), bd)
      p3, p4 = joinpath(dir, "mps_stop.txt"), joinpath(dir, "mps_end.txt")
      mps_circuit_dynamics(basis_mps(N, ones), gates, phases, 1e-15, O, p3; layer_ends = rec, χmax = χm)
      mps_circuit_dynamics(basis_mps(N, ones), gates, phases, 1e-15, O, p4; layer_ends = rec)
      @test lm < length(rec)                    # it does stop early
      @test length(stoplines(p3)) == 1 && occursin("stopped at gate $(rec[lm]) ", only(stoplines(p3))) &&
            occursin("bond dim", only(stoplines(p3)))
      @test length(stoplines(p4)) == 1 && occursin("stopped at gate $M ", only(stoplines(p4))) &&
            occursin("end of circuit", only(stoplines(p4)))

      # showprogress = true only adds progress output: the same stop and values
      quiet(f) = redirect_stdio(f; stdout = devnull, stderr = devnull)
      _, ss, evss, = campssrc_circuit_dynamics(ψ0, gates, phases, 1e-11, O; layer_ends = rec, χmax = χ, track = true)
      _, sp, evsp, = quiet(() -> campssrc_circuit_dynamics(ψ0, gates, phases, 1e-11, O;
        layer_ends = rec, χmax = χ, track = true, showprogress = true))
      @test ss == sp == rec[l] && evsp == evss
      _, evsm = mps_circuit_dynamics(basis_mps(N, ones), gates, phases, 1e-15, O, nothing; layer_ends = rec, χmax = χm)
      _, evsmp = quiet(() -> mps_circuit_dynamics(basis_mps(N, ones), gates, phases, 1e-15, O, nothing;
        layer_ends = rec, χmax = χm, showprogress = true))
      @test length(evsm) == lm + 1 && evsmp == evsm
    end
  end

  @testset "χmax: the configurations tell the rules apart ($k)" for k in sort(collect(keys(seen)))
    @test seen[k]
  end
end
