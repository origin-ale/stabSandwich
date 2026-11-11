# TimeEvolution.jl against the dense reference at small N, and against exact classical
# results at N = 46. Expected values: CONVENTIONS.md and DenseReference.jl only.
#
# Per-method parameters are kept separate on purpose (χ_camps, χ_campspp, thl_pp, …).

using ITensors, ITensorMPS

# -- Test circuits ----------------------------------------------------------------------

"Doped π/4 XXZ circuit as in the production runs: `:xy` doping at 3π/16 on a seeded mask."
function doped_xxz(N, t; p = 0.3, seed = 11, magicphase = 3π/16, kind = :xy)
  gates, phases = xxz_circuit(π/4, π/4, t, N)
  mask = doping_mask(Xoshiro(seed), phases, p)
  phases = kind == :xy ? xy_magic(phases, mask; magicphase) : z_magic(phases, mask; magicphase)
  return gates, phases, mask
end

basis_mps(N, ones) = MPS(ComplexF64, siteinds("Qubit", N), [j in ones ? "1" : "0" for j in 1:N])

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

# The three small configurations every method is checked on.
function small_configs()
  cfgs = []
  # A: doped π/4 XXZ (Clifford π/4 gates plus general 3π/16 gates), layer by layer.
  N, t, ones = 6, 4, [2, 4, 5, 6]
  gates, phases, mask = doped_xxz(N, t)
  @assert 2 <= count(mask) < length(mask)
  push!(cfgs, (name = "doped π/4 XXZ, N=6", N, ones, gates, phases, layer_ends = spec_layer_ends(N, t)))
  # B: XXZ with general angles everywhere, odd N.
  N, t, ones = 5, 3, [1, 4]
  gates, phases = xxz_circuit(0.3, 0.7, t, N)
  push!(cfgs, (name = "general XXZ, N=5", N, ones, gates, phases, layer_ends = spec_layer_ends(N, t)))
  # C: random high-weight rotations, recorded after every gate (no layer ends).
  N, ones = 4, [3]
  gates, phases = rotation_circuit(Xoshiro(21), 10, N)
  @assert all(abs(ϕ - c) > 1e-3 for ϕ in phases, c in (0:8) .* (π / 4))
  push!(cfgs, (name = "random rotations, N=4", N, ones, gates, phases, layer_ends = nothing))
  return cfgs
end

record_points(cfg) = isnothing(cfg.layer_ends) ? collect(eachindex(cfg.gates)) : cfg.layer_ends

@testset "TimeEvolution" begin
  switched = Bool[]      # did CAMPS-PP hand over to PP? (collected over all configurations)

  @testset "exact regime vs dense: $(cfg.name)" for cfg in small_configs()
    (; N, ones, gates, phases, layer_ends) = cfg
    M = length(gates)
    O = transferredmagnetization(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, record_points(cfg))
    @test ref[1] ≈ 0 atol = EXACT_ATOL
    @test maximum(abs, ref) > 0.1                 # the dynamics is not trivial
    ψ0 = computationalcamps(N, ones)

    @testset "camps_circuit_dynamics" begin
      χ_camps = 10^6                              # never reached: runs the whole circuit
      ψ, k, s, evs = camps_circuit_dynamics(ψ0, gates, phases, χ_camps, O; layer_ends)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test 0 <= k <= N
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
      @test 0 <= s_switch <= M
      @test 0 <= k_switch <= N
      # Frame-dependent: whether and where the switch happens.
      χ_campspp == 10^6 ? (@test s_switch == M) : push!(switched, s_switch < M)
    end

    @testset "campssrc_circuit_dynamics ($crit, $strat)" for (crit, strat) in
        [(:chi3, :snake), (:entangle, :brickwork)]
      thl_camps = 1e-15                           # SVD cutoff: nothing is discarded
      ψ, s, evs = campssrc_circuit_dynamics(ψ0, gates, phases, thl_camps, O;
        criterion = crit, strategy = strat, layer_ends)
      @test s == M
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
    end

    @testset "mps_circuit_dynamics" begin
      thl_mps = 1e-15
      ψ, evs = mps_circuit_dynamics(basis_mps(N, ones), gates, phases, thl_mps, O, nothing; layer_ends)
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
    end
  end

  @testset "CAMPS-PP switches to PP in the small configurations" begin
    # Guard against a vacuous pass: the χ = 2 and χ = 4 runs above must include real
    # hand-overs (a configuration whose bond dimension never reaches χ does not count).
    @test count(switched) >= 2
  end

  @testset "CAMPS-PP continues from a CAMPS state" begin
    # Evolve s gates with CAMPS, hand over to PP for the rest: the result must be the
    # dense value at the end, wherever the hand-over is.
    # N = 4 keeps the untruncated Pauli sums small (at most 4^4 terms).
    N, t, ones = 4, 3, [3, 4]
    gates, phases, _ = doped_xxz(N, t; p = 0.5)
    M = length(gates)
    O = transferredmagnetization(N, ones)
    ref = dense_dynamics(N, ones, gates, phases, O, [M])
    for s0 in (0, 5, 14, M - 1, M)
      ψ, k = evolve(computationalcamps(N, ones), s0, gates, phases)
      s, evs = pauliprop_circuit_dynamics(ψ, s0, gates, phases, 0.0, typemax(Int), O;
        layer_ends = [M])
      @test s == M
      if s0 < M
        @test length(evs) == 1
        @test evs[1] ≈ ref[2] atol = EXACT_ATOL
      else
        @test isempty(evs)
      end
    end
  end

  @testset "truncated regime: PP threshold" begin
    # A truncated run must approach the dense result as the threshold is tightened.
    # No absolute bound is asserted for the loose thresholds (none is stated in the spec).
    N, t, ones = 6, 4, [2, 4, 5, 6]
    gates, phases, _ = doped_xxz(N, t; p = 0.5)
    O = transferredmagnetization(N, ones)
    le = spec_layer_ends(N, t)
    ref = dense_dynamics(N, ones, gates, phases, O, le)
    errs = map((3e-1, 3e-2, 1e-3, 1e-6, 0.0)) do thl_pp
      _, evs = pauliprop_circuit_dynamics(ones, gates, phases, thl_pp, typemax(Int), O; layer_ends = le)
      maxabsdiff(evs, ref)
    end
    @test errs[end] <= EXACT_ATOL
    @test errs[1] > 1e-3                          # the loosest threshold really truncates
    @test issorted(collect(errs); rev = true)     # error does not grow as thl shrinks
  end

  @testset "truncated regime: Nmax stops PP early" begin
    N, t, ones = 6, 4, [2, 4, 5, 6]
    gates, phases, _ = doped_xxz(N, t; p = 0.5)
    O = transferredmagnetization(N, ones)
    le = spec_layer_ends(N, t)
    ref = dense_dynamics(N, ones, gates, phases, O, le)
    Nmax_pp = 20
    s, evs, nterms = pauliprop_circuit_dynamics(ones, gates, phases, 0.0, Nmax_pp, O;
      layer_ends = le, track = true)
    @test s < length(gates)                       # stopped before the end
    @test s in le                                 # … at a layer end
    @test nterms[end] >= Nmax_pp && all(nterms[1:end-1] .< Nmax_pp)
    @test maxabsdiff(evs, ref[1:length(evs)]) <= EXACT_ATOL   # what was computed is exact
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
