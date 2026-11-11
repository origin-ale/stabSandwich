# TimeEvolution.jl against the dense reference at small N, and against exact classical
# results at N = 46. Expected values: CONVENTIONS.md and DenseReference.jl only.
#
# Per-method parameters are kept separate on purpose (χ_camps, χ_campspp, thl_pp, …).

using ITensors, ITensorMPS

const bonddim = DisentangleCAMPS.bonddim

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

  late_switch = Bool[]   # CAMPS-PP hand-overs that happen after the first record point
  early_visible = Bool[] # would a value recorded two gates early differ from the reference?
  late_visible = Bool[]  # … and one gate late?
  # Did the layer-end disentangler change the frame? Per (criterion, strategy).
  frame_changed = Dict((:chi3, :snake) => Bool[], (:entangle, :brickwork) => Bool[])

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
    ψ0 = computationalcamps(N, ones)

    @testset "camps_circuit_dynamics" begin
      χ_camps = 10^6                              # never reached: runs the whole circuit
      ψ, k, s, evs = camps_circuit_dynamics(ψ0, gates, phases, χ_camps, O; layer_ends)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test 0 <= k <= N
      ev = cmps.expectation(ψ, O)                 # CONVENTIONS.md §5: a real number
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

    @testset "campspp_circuit_dynamics, χ = $χ_campspp" for χ_campspp in (2, 4, 10^6)
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
          @test bonddim(evolve(computationalcamps(N, ones), s_switch - 1, gates, phases)[1]) < χ_campspp
          push!(late_switch, s_switch > first(rec))
        end
      end
    end

    @testset "campssrc_circuit_dynamics ($crit, $strat)" for (crit, strat) in
        [(:chi3, :snake), (:entangle, :brickwork)]
      thl_camps = 1e-15                           # SVD cutoff: nothing is discarded
      ψ, s, evs = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
        criterion = crit, strategy = strat, layer_ends)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test abs(imag(cmps.expectation(ψ, O))) <= EXACT_ATOL
      # Frame-dependent: compare the frame with the same gates applied without the
      # layer-end disentangler.
      plain = computationalcamps(N, ones)
      for (g, ϕ) in zip(gates, phases)
        apply!(plain, PauliOperator(getpauli(g, N)), ϕ, thl_camps)
      end
      push!(frame_changed[(crit, strat)], ψ.Cdag != plain.Cdag)
    end

    @testset "mps_circuit_dynamics" begin
      thl_mps = 1e-15
      ψ, evs = mps_circuit_dynamics(basis_mps(N, ones), gates, phases, thl_mps, O, nothing; layer_ends)
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
    end
  end

  @testset "the small configurations exercise the hand-over and the disentangler" begin
    # Guards against vacuous passes above (frame-dependent facts about the old code).
    @test count(late_switch) >= 2       # hand-overs after the first record point
    @test count(early_visible) >= 2     # record points two gates early would be noticed
    @test count(late_visible) >= 2      # … and one gate late
    # The layer-end disentangler does change the frame with (:entangle, :brickwork).
    @test count(frame_changed[(:entangle, :brickwork)]) >= 2
    # With (:chi3, :snake), the production setting (called with min_diff = 6), it never
    # acts at these sizes: the exact-regime checks above cover that path only as a no-op.
    # Its real coverage is the N = 46 old-vs-new replay (plan, T11). If this starts to
    # fail, the small configurations have begun to exercise it: turn it into a positive check.
    @test !any(frame_changed[(:chi3, :snake)])
  end

  @testset "CAMPS-PP continues from a CAMPS state, layer by layer" begin
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
      ψ, k = evolve(computationalcamps(N, ones), s0, gates, phases)
      s, evs = pauliprop_circuit_dynamics(ψ, s0, gates, phases, 0.0, typemax(Int), O; layer_ends = le)
      later = [l for l in 1:t if le[l] > s0]                      # layers still to come
      @test s == M
      @test length(evs) == length(later)
      isempty(later) || @test maxabsdiff(evs, ref[later .+ 1]) <= EXACT_ATOL
    end
  end

  @testset "CAMPS-PP: `thl` and `Nmax` reach the PP stage" begin
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    M = length(gates)
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    ψ0 = computationalcamps(N, ones)
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

  @testset "truncated regime: CAMPS SVD threshold" begin
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    ψ0 = computationalcamps(N, ones)
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

  @testset "output file" begin
    # With a file name as `output`, every recorded value is appended as "index<TAB>value".
    cfg = small_configs()[1]
    (; N, ones, gates, phases, layer_ends) = cfg
    O = test_observable(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, layer_ends)
    ψ0 = computationalcamps(N, ones)
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
      ψ0 = computationalcamps(N, ones)
      O = transferredmagnetization(N, ones)
      le = spec_layer_ends(N, t)
      exact = swap_dynamics(N, ones, t)
      gates, undoped = xxz_circuit(π/4, π/4, t, N)
      M = length(gates)
      χ_campspp, thl_campspp, Nmax_campspp = 64, 1e-10, 10^6
      thl_camps = 1e-12

      run_pp(phases) = campspp_circuit_dynamics(ψ0, χ_campspp, thl_campspp, Nmax_campspp,
        gates, phases, O; layer_ends = le, track = true)
      run_src(phases) = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
        criterion = :chi3, strategy = :snake, layer_ends = le, track = true)

      evs0, _, bds0, _, _, _, s_switch0, k0 = run_pp(undoped)
      @test maxabsdiff(evs0, exact) <= EXACT_ATOL  # undoped run against the exact SWAP result
      _, _, src0, srcbds0 = run_src(undoped)
      @test maxabsdiff(src0, exact) <= EXACT_ATOL

      for α in (π/3, 1.234)
        _, phases, mask = doped_xxz(N, t; p = 0.3, kind = :z, magicphase = α)
        @assert count(mask) > 10
        evs, s, bds, _, _, _, s_switch, k = run_pp(phases)
        @test evs == evs0                          # bitwise equal to the p = 0 run
        @test s == M && s_switch == M              # no switch to PP
        @test all(==(1), bds) && k == 0            # frame-dependent: nothing was added
        # The truncating method sends every doped gate through an SVD and never
        # renormalises, and `expectation` is not normalised (CONVENTIONS.md §5). Over a few
        # hundred doped gates the norm drifts by rounding (measured 4e-12 on norm² at 23
        # layers), which scales ⟨O⟩ ≈ 30 by the same factor. So: the norm stays within 1e-10
        # of 1, each layer agrees up to that drift, and the normalised final value agrees.
        ψsrc, _, src, srcbds = run_src(phases)
        drift = abs(norm(ψsrc)^2 - 1)
        @test drift <= EXACT_ATOL
        @test all(abs.(src .- src0) .<= EXACT_ATOL .+ abs.(src0) .* drift)
        @test abs(src[end] / norm(ψsrc)^2 - src0[end]) <= EXACT_ATOL
        @test all(==(1), srcbds)                   # frame-dependent
      end
    end
  end

  @testset "large N, Clifford only: undoped π/4 XXZ at N=46" begin
    N, t = 46, get(ENV, "STABSANDWICH_SLOW", "") == "" ? 6 : 23
    _, ones = domainwallstate(Xoshiro(100), N, 0.6)
    O = transferredmagnetization(N, ones)
    le = spec_layer_ends(N, t)
    exact = swap_dynamics(N, ones, t)
    @test maximum(abs, exact) >= 2                # excitations do cross the cut
    gates, phases = xxz_circuit(π/4, π/4, t, N)
    ψ0 = computationalcamps(N, ones)

    ψ, k, s, evs = camps_circuit_dynamics(ψ0, gates, phases, 64, O; layer_ends = le)
    @test s == length(gates) && k == 0
    @test maxabsdiff(evs, exact) <= EXACT_ATOL
    @test DisentangleCAMPS.bonddim(ψ) == 1         # frame-dependent

    s, evs = pauliprop_circuit_dynamics(ones, gates, phases, 1e-10, 10^6, O; layer_ends = le)
    @test s == length(gates)
    @test maxabsdiff(evs, exact) <= EXACT_ATOL
  end
end
