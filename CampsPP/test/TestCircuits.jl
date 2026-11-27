# Circuits.jl against CONVENTIONS.md §2. Convention-invariant (no frame involved).

@testset "Circuits" begin
  @testset "xxz_circuit structure, N=$N" for N in (2, 4, 5, 6)
    t = 3
    ϕxy, ϕz = 0.3, 0.7
    rots, phases = xxz_circuit(ϕxy, ϕz, t, N)
    bonds = spec_bonds(N)
    @test length(bonds) == N - 1
    @test length(rots) == 3 * (N - 1) * t
    @test length(phases) == length(rots)

    ok_syms = true; ok_inds = true; ok_phases = true
    for l in 1:t, (b, bond) in enumerate(bonds)
      m = (l - 1) * (N - 1) + b                   # m-th triple in gate order
      for (g, (syms, ϕ)) in enumerate([([:X, :X], ϕxy), ([:Y, :Y], ϕxy), ([:Z, :Z], ϕz)])
        i = 3m - 3 + g
        ok_syms &= collect(rots[i].symbols) == syms
        ok_inds &= collect(rots[i].qinds) == bond
        ok_phases &= phases[i] == ϕ
      end
    end
    @test ok_syms
    @test ok_inds
    @test ok_phases

    @test layerends(N, t, xxz_circuit) == spec_layer_ends(N, t)
  end

  @testset "layer of 135 gates at N=46" begin
    @test layerends(46, 23, xxz_circuit) == collect(135:135:3105)
    @test length(xxz_circuit(π/4, π/4, 23, 46)[1]) == 3105
  end

  @testset "XXZ triple at ϕ_xy = ϕ_z = π/4 is a SWAP (dense)" begin
    # exp(iπ/4 (XX+YY+ZZ)) = e^{iπ/4} SWAP, so basis states get their two bits exchanged.
    N = 4
    rots, phases = xxz_circuit(π/4, π/4, 1, N)
    for ones in ([1], [2, 3], [1, 2, 4], Int[])
      bits = [j in ones for j in 1:N]
      ψ = DR.basis_state(N, ones)
      for (m, bond) in enumerate(spec_bonds(N))
        for i in 3m-2:3m
          DR.rotate!(ψ, gatestring(rots[i], N), phases[i])
        end
        bits[bond[1]], bits[bond[2]] = bits[bond[2]], bits[bond[1]]
        @test all(DR.bloch(ψ, j)[3] ≈ (bits[j] ? -1.0 : 1.0) for j in 1:N)
      end
    end
  end

  @testset "doping with a mask" begin
    N, t = 5, 2
    _, phases = xxz_circuit(π/4, π/4, t, N)
    ntriples = (N - 1) * t
    mask = falses(ntriples); mask[[2, 5, 8]] .= true
    α = 3π/16

    xy = xy_magic(phases, mask; magicphase = α)
    z = z_magic(phases, mask; magicphase = α)
    @test length(xy) == length(phases) == length(z)
    @test all(phases .== π/4)                      # the input is not modified
    for m in 1:ntriples
      if mask[m]
        @test xy[3m-2] == α && xy[3m-1] == α && xy[3m] == π/4
        @test z[3m-2] == π/4 && z[3m-1] == π/4 && z[3m] == α
      else
        @test xy[3m-2] == π/4 && xy[3m-1] == π/4 && xy[3m] == π/4
        @test z[3m-2] == π/4 && z[3m-1] == π/4 && z[3m] == π/4
      end
    end

    @test xy_magic(phases, falses(ntriples); magicphase = α) == phases
    @test all(z_magic(phases, trues(ntriples); magicphase = α)[3:3:end] .== α)
  end

  @testset "doping_mask" begin
    N, t = 46, 23
    _, phases = xxz_circuit(π/4, π/4, t, N)
    n = (N - 1) * t                                # one entry per triple, whole circuit
    @test n == 1035

    mask = doping_mask(Xoshiro(100), phases, 0.04)
    @test length(mask) == n
    @test eltype(mask) == Bool
    @test mask == doping_mask(Xoshiro(100), phases, 0.04)      # reproducible from the seed
    @test mask != doping_mask(Xoshiro(101), phases, 0.04)
    @test !any(doping_mask(Xoshiro(1), phases, 0.0))
    @test all(doping_mask(Xoshiro(1), phases, 1.0))

    # Each entry is true with probability p: the count is Binomial(n, p). Fixed seeds, 5σ.
    for p in (0.01, 0.04, 0.3)
      counts = [count(doping_mask(Xoshiro(seed), phases, p)) for seed in 1:40]
      σ = sqrt(n * p * (1 - p))
      @test all(abs.(counts .- n * p) .<= 5σ)
      @test abs(sum(counts) / 40 - n * p) <= 5σ / sqrt(40)
    end
  end

  @testset "doping_mask: every position is doped with probability p" begin
    _, phases = xxz_circuit(π/4, π/4, 4, 12)
    n = length(phases) ÷ 3
    p, seeds = 0.3, 400
    freq = zeros(n)
    for seed in 1:seeds
      freq .+= doping_mask(Xoshiro(seed), phases, p) ./ seeds
    end
    σ = sqrt(p * (1 - p) / seeds)
    @test all(abs.(freq .- p) .<= 5σ)              # no position is favoured
  end

  @testset "doping_mask: one rand(rng) per triple, in order (§2 probability mode)" begin
    # The T11 replay regenerates the Phase 0 masks with MersenneTwister(hash((seed + s, p))).
    # Pinned for MersenneTwister only, the RNG the scripts and the replay use.
    _, phases = xxz_circuit(π/4, π/4, 23, 46)
    T = length(phases) ÷ 3
    for p in (0.01, 0.04, 0.5), s in 0:3
      seed = hash((100 + s, p))
      r = MersenneTwister(seed)
      draws = [rand(r) for _ in 1:T]
      @test doping_mask(MersenneTwister(seed), phases, p) == (draws .< p)
    end
  end

  # -- Post-port item 9: count mode (CONVENTIONS.md §2 "Doping") --------------------------
  # Expected values from the spec: exactly n distinct triples, drawn uniformly without
  # replacement, so each triple is doped with probability n/T, each pair with
  # n(n−1)/(T(T−1)), and every n-subset is equally likely. Frequency checks are 5σ over
  # fixed seeds; each helper has a negative check against a wrong sampler.
  if !isdefined(CampsPP, :doping_mask_count)
    @testset "doping_mask_count [item 9, not implemented yet]" begin
      @test_broken isdefined(CampsPP, :doping_mask_count)
    end
  else
    # Masks from `sampler(seed)` for seeds 1:S.
    masks_over(sampler, S) = [sampler(seed) for seed in 1:S]

    # Is every position's frequency within 5σ of q = n/T?
    function positions_uniform(masks, q)
      S = length(masks)
      freq = sum(masks) ./ S
      return all(abs.(freq .- q) .<= 5 * sqrt(q * (1 - q) / S))
    end

    # Is every pair's co-doping frequency within 5σ of q2?
    function pairs_uniform(masks, q2)
      S, T = length(masks), length(first(masks))
      M = Int.(reduce(hcat, masks))                # T × S
      C = (M * M') ./ S
      σ = sqrt(q2 * (1 - q2) / S)
      return all(abs(C[a, b] - q2) <= 5σ for a in 1:T for b in a+1:T)
    end

    # Is every n-subset of 1:T drawn, each with frequency within 5σ of 1/binomial(T, n)?
    function subsets_uniform(masks, T, n)
      S, K = length(masks), binomial(T, n)
      counts = Dict{Vector{Int}, Int}()
      for m in masks
        k = findall(m)
        counts[k] = get(counts, k, 0) + 1
      end
      q = 1 / K
      σ = sqrt(q * (1 - q) / S)
      return length(counts) == K && all(abs(c / S - q) <= 5σ for c in values(counts))
    end

    # Wrong samplers for the negative checks (each returns a T-entry Bool mask)
    first_n(T, n) = seed -> (m = falses(T); m[1:n] .= true; m)
    cyclic_block(T, n) = seed -> begin                # uniform marginals, exact count
      s = rand(Xoshiro(seed), 0:T-1); m = falses(T); m[mod1.(s .+ (1:n), T)] .= true; m
    end
    iid(T, n) = seed -> rand(Xoshiro(seed), T) .< n / T   # independent draws, p = n/T
    never_last(T, n) = seed -> (m = falses(T); m[randperm(Xoshiro(seed), T - 1)[1:n]] .= true; m)

    count_mask(phases, n) = seed -> doping_mask_count(Xoshiro(seed), phases, n)

    @testset "doping_mask_count: exactly n distinct triples" begin
      for (N, t) in ((46, 23), (5, 2), (2, 1))
        _, phases = xxz_circuit(π/4, π/4, t, N)
        T = (N - 1) * t
        ns = unique(filter(n -> 0 <= n <= T, [0, 1, 2, T ÷ 2, 41, T - 1, T]))
        ok = true
        for n in ns, seed in 1:5, rng in (Xoshiro(seed), MersenneTwister(seed))
          m = doping_mask_count(rng, phases, n)
          ok &= length(m) == T && eltype(m) == Bool && count(m) == n
        end
        @test ok
        @test !any(doping_mask_count(Xoshiro(3), phases, 0))
        @test all(doping_mask_count(Xoshiro(3), phases, T))
      end
    end

    @testset "doping_mask_count: invalid n" begin
      _, phases = xxz_circuit(π/4, π/4, 2, 5)         # T = 8
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, -1)
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, 9)
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, 2.5)
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, typemax(Int))
      @test count(doping_mask_count(Xoshiro(1), phases, 8)) == 8   # the bound itself is valid
      # §2: n must be a Julia `Integer`; floats with an integral value are refused (review
      # follow-up 3; mutant: accept any n with isinteger(n)). 3//1 is integral too.
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, 3.0)
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, 8.0)
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, 3//1)
    end

    @testset "doping_mask_count: any Integer type, same mask as Int" begin
      # §2: n is any Julia `Integer` (not a `Bool`), and the mask depends only on the RNG
      # state, so a non-Int n gives the mask of the equal Int (review follow-up 1; mutant:
      # only `Int` accepted).
      _, phases = xxz_circuit(π/4, π/4, 23, 46)        # T = 1035
      ok = true
      for n in (0, 1, 41, 1035), conv in (Int32, UInt, BigInt), seed in 1:3
        for mk in (Xoshiro, MersenneTwister)
          m = doping_mask_count(mk(seed), phases, conv(n))
          ok &= count(m) == n && m == doping_mask_count(mk(seed), phases, n)
        end
      end
      @test ok
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, UInt(1036))
      @test_throws ArgumentError doping_mask_count(Xoshiro(1), phases, big(-1))
    end

    @testset "doping_mask_count: advances rng; the 2-argument form uses the global RNG" begin
      # §2: the call draws from `rng` itself and leaves it advanced, so successive calls on
      # one rng give independent masks (review follow-up 2; mutant: draw from copy(rng)).
      # Two independent 41-subsets of 1035 are equal with probability 1/binomial(1035, 41),
      # below 1e-70, so "differ" is not luck.
      _, phases = xxz_circuit(π/4, π/4, 23, 46)
      n = 41
      for mk in (Xoshiro, MersenneTwister)
        r = mk(4); r0 = copy(r)
        a = doping_mask_count(r, phases, n)
        @test r != r0                                  # the state moved
        b = doping_mask_count(r, phases, n)
        @test a != b                                   # a fresh mask, not the first again
        @test a == doping_mask_count(mk(4), phases, n) # the first call is the seeded mask
      end
      # Docstring: `doping_mask_count(phases, n)` draws from the global RNG (mutant: a fixed
      # Xoshiro(0)). After the same seed it equals the 3-argument call on default_rng(),
      # successive calls differ, and another seed gives another mask.
      Random.seed!(11); a = doping_mask_count(phases, n)
      Random.seed!(11); b = doping_mask_count(Random.default_rng(), phases, n)
      @test count(a) == n && a == b
      Random.seed!(11); c1 = doping_mask_count(phases, n); c2 = doping_mask_count(phases, n)
      @test c1 == a && c1 != c2
      Random.seed!(12); @test doping_mask_count(phases, n) != a
    end

    @testset "doping_mask_count: depends only on the RNG state" begin
      _, phases = xxz_circuit(π/4, π/4, 23, 46)
      T, n = length(phases) ÷ 3, 41
      m = doping_mask_count(MersenneTwister(hash((101, n))), phases, n)
      @test m == doping_mask_count(MersenneTwister(hash((101, n))), phases, n)
      @test m != doping_mask_count(MersenneTwister(hash((102, n))), phases, n)
      @test m != doping_mask_count(MersenneTwister(hash((101, n + 1))), phases, n)
      # Same state reached another way (a copy), and phases with other values: same mask
      r = MersenneTwister(7); rand(r, 13)
      @test doping_mask_count(copy(r), phases, n) == doping_mask_count(r, phases, n)
      @test doping_mask_count(Xoshiro(5), rand(Xoshiro(9), 3T), n) ==
            doping_mask_count(Xoshiro(5), phases, n)
    end

    @testset "doping_mask_count: every triple doped with probability n/T" begin
      _, phases = xxz_circuit(π/4, π/4, 4, 12)        # T = 44
      T, S = length(phases) ÷ 3, 4000
      for n in (1, 11, 30)
        @test positions_uniform(masks_over(count_mask(phases, n), S), n / T)
      end
      @test !positions_uniform(masks_over(first_n(T, 11), S), 11 / T)          # negative
      @test !positions_uniform(masks_over(never_last(T, 11), S), 11 / T)       # negative
    end

    @testset "doping_mask_count: first and last triple unbiased (N=46)" begin
      _, phases = xxz_circuit(π/4, π/4, 23, 46)
      T, n, S = length(phases) ÷ 3, 41, 4000
      masks = masks_over(count_mask(phases, n), S)
      q, σ = n / T, sqrt(n / T * (1 - n / T) / S)
      firstlast_ok(ms) = all(abs(sum(m[k] for m in ms) / S - q) <= 5σ for k in (1, T))
      @test firstlast_ok(masks)
      @test !firstlast_ok(masks_over(first_n(T, n), S))                         # negative
      @test !firstlast_ok(masks_over(never_last(T, n), S))                      # negative
    end

    @testset "doping_mask_count: pairs at n(n−1)/(T(T−1)), not (n/T)²" begin
      _, phases = xxz_circuit(π/4, π/4, 2, 6)         # T = 10
      T, n, S = length(phases) ÷ 3, 5, 20000
      q2 = n * (n - 1) / (T * (T - 1))                # 2/9, vs (n/T)² = 1/4
      @test pairs_uniform(masks_over(count_mask(phases, n), S), q2)
      @test !pairs_uniform(masks_over(iid(T, n), S), q2)                        # negative
      @test !pairs_uniform(masks_over(cyclic_block(T, n), S), q2)               # negative
      @test positions_uniform(masks_over(cyclic_block(T, n), S), n / T)         # (marginals fine)
    end

    @testset "doping_mask_count: every n-subset equally likely" begin
      _, phases = xxz_circuit(π/4, π/4, 2, 4)         # T = 6
      T, S = length(phases) ÷ 3, 20000
      for n in (2, 3)
        @test subsets_uniform(masks_over(count_mask(phases, n), S), T, n)
      end
      @test !subsets_uniform(masks_over(cyclic_block(T, 3), S), T, 3)           # negative
    end

    @testset "count mask applied at N=4 (end to end)" begin
      # Hand pin. N=4, t=2: bonds (1,2), (3,4), (2,3) per layer, so triples 1…6 sit on
      # bonds (1,2), (3,4), (2,3), (1,2), (3,4), (2,3). Triples {2, 6}: XX at gates 4, 16,
      # YY at 5, 17 (:xy); ZZ at 6, 18 (:z).
      rots, phases = xxz_circuit(π/4, π/4, 2, 4)
      α = 3π/16
      mask = Bool[0, 1, 0, 0, 0, 1]
      @test findall(xy_magic(phases, mask; magicphase = α) .!= phases) == [4, 5, 16, 17]
      @test findall(z_magic(phases, mask; magicphase = α) .!= phases) == [6, 18]
      @test [gatestring(rots[i], 4) for i in (4, 5, 6, 16, 17, 18)] ==
            ["IIXX", "IIYY", "IIZZ", "IXXI", "IYYI", "IZZI"]

      # Every count mask: exactly n XX-YY pairs (2n gates) or n ZZ gates, at 3m−2, 3m−1 / 3m
      ok = true
      for n in 0:6, seed in 1:10
        m = doping_mask_count(MersenneTwister(seed), phases, n)
        ks = findall(m)
        xy = findall(xy_magic(phases, m; magicphase = α) .!= phases)
        z = findall(z_magic(phases, m; magicphase = α) .!= phases)
        ok &= length(ks) == n && xy == sort([3ks .- 2; 3ks .- 1]) && z == 3ks
        ok &= all(gatestring(rots[i], 4) in ("XXII", "IIXX", "IXXI") for i in 3ks .- 2)
        ok &= all(gatestring(rots[i], 4) in ("YYII", "IIYY", "IYYI") for i in 3ks .- 1)
      end
      @test ok
    end
  end

  @testset "rotation_circuit" begin
    N, t = 5, 12
    rots, phases = rotation_circuit(Xoshiro(7), t, N)
    @test length(rots) == t == length(phases)
    @test all(collect(r.qinds) == collect(1:N) for r in rots)   # every gate acts on all N qubits
    @test all(0 .<= phases .< 2π)
    rots2, phases2 = rotation_circuit(Xoshiro(7), t, N)
    @test [gatestring(r, N) for r in rots] == [gatestring(r, N) for r in rots2]
    @test phases == phases2

    # Given strings and angles are returned as they are.
    Ps = [[:X, :I, :Z], [:Y, :Y, :I]]
    ϕs = [0.2, -1.3]
    rots, phases = rotation_circuit(Ps, ϕs, 3)
    @test [gatestring(r, 3) for r in rots] == ["XIZ", "YYI"]
    @test phases == ϕs
  end

  # -- Kicked Ising (CONVENTIONS.md §2 "Kicked Ising"; KIM plan step 2) -------------------
  # Expected values from the spec (KimSpec.jl builds the circuit from it), hand derivations
  # and the dense reference. The spec circuit's own checks run alongside those of
  # `kicked_ising_circuit` and `kicked_ising_layerends`.

  # The builder's circuit as (strings, phases), for the dense reference
  function kim_builder_circuit(ε, t, N)
    local rots, phases = kicked_ising_circuit(ε, t, N)
    return [gatestring(r, N) for r in rots], phases
  end
  kim_sources = Any[("spec circuit", kim_spec_circuit), ("kicked_ising_circuit", kim_builder_circuit)]

  @testset "kicked_ising_circuit structure, N=$N" for N in (1, 2, 3, 5, 6)
    t = 3
    L = kim_spec_L(N)
    @test L == Dict(1 => 2, 2 => 5, 3 => 9, 5 => 15, 6 => 18)[N]    # §2: 2N + B
    bonds = kim_spec_bonds(N)
    @test length(bonds) == (N >= 3 ? N : N - 1)
    # Distinct ε_j on every qubit, so ε taken from the wrong qubit shows
    εv = [0.1j - 0.037 for j in 1:N]
    rots, phases = kicked_ising_circuit(εv, t, N)
    @test length(rots) == L * t
    @test length(phases) == L * t
    @test phases isa AbstractVector{Float64}

    # Gate i of layer l at position (l − 1)·L + i: ZZ on the bonds (bricklayer order, then
    # [1, N]), then Z_1…Z_N, then Y_1…Y_N; angles −π/4, −(π/2 + ε_j)/2, −π/4.
    want = [[(:ZZ, b, -π/4) for b in bonds];
            [(:Z, [j], -(π/2 + εv[j]) / 2) for j in 1:N];
            [(:Y, [j], -π/4) for j in 1:N]]
    @test length(want) == L
    ok_syms = true; ok_inds = true; ok_phases = true
    for l in 1:t, (i, (sym, qs, ϕ)) in enumerate(want)
      g = (l - 1) * L + i
      ok_syms &= collect(rots[g].symbols) == (sym == :ZZ ? [:Z, :Z] : [sym])
      ok_inds &= collect(rots[g].qinds) == qs
      ok_phases &= isapprox(phases[g], ϕ; atol = 1e-14)
    end
    @test ok_syms
    @test ok_inds
    @test ok_phases
    # The wrap bond is the last ZZ gate of each layer, listed as [1, N], for N ≥ 3 only;
    # N = 2 has exactly one ZZ gate per layer, N = 1 none
    @test count(r -> collect(r.symbols) == [:Z, :Z], rots[1:L]) == (N >= 3 ? N : N - 1)
    N >= 3 && @test collect(rots[N].qinds) == [1, N]
    # The non-Clifford angles are exactly the Z ones (one per qubit per layer)
    @test count(ϕ -> !isapprox(ϕ, -π/4; atol = 1e-12), phases) == N * t

    # The number form equals the vector form with fill(ε, N)
    for ε in (0.0, 0.3, -0.11)
      r1, p1 = kicked_ising_circuit(ε, t, N)
      r2, p2 = kicked_ising_circuit(fill(ε, N), t, N)
      @test [gatestring(r, N) for r in r1] == [gatestring(r, N) for r in r2]
      @test [collect(r.qinds) for r in r1] == [collect(r.qinds) for r in r2]
      @test p1 == p2
    end
    # The builder's circuit is the spec circuit, gate by gate
    strings, ph = kim_builder_circuit(εv, t, N)
    spec_strings, spec_ph = kim_spec_circuit(εv, t, N)
    @test strings == spec_strings
    @test maxabsdiff(ph, spec_ph) <= 1e-14

    # Any other vector length is an ArgumentError
    @test_throws ArgumentError kicked_ising_circuit(fill(0.1, N - 1), t, N)
    @test_throws ArgumentError kicked_ising_circuit(fill(0.1, N + 1), t, N)

    # Layer ends L·l, L from the spec formula (not from the builder)
    @test kicked_ising_layerends(N, t) == [L * l for l in 1:t]
  end

  @testset "kicked_ising_circuit at N=46: 138 gates per layer" begin
    @test kicked_ising_layerends(46, 23) == collect(138:138:3174)
    rots, phases = kicked_ising_circuit(0.0, 23, 46)
    @test length(rots) == length(phases) == 3174
    @test collect(rots[46].qinds) == [1, 46]           # the wrap bond closes the ZZ block
    @test collect(rots[47].symbols) == [:Z] && collect(rots[138].symbols) == [:Y]
  end

  # Measured (dense, 2026-11-26; N = 4…14, checked by an independent review): from the Néel
  # state, ⟨Z_j⟩ = 0 after layers 1 … ⌈N/2⌉ for every ε ≠ 0 (also per-site ε, also from
  # |0…0⟩); it first differs from 0 at layer ⌈N/2⌉ + 1, once the light cone has wrapped
  # around the ring (on the open chain: one layer later for even N, the same layer for odd
  # N). At ε = 0 it can stay 0 much longer (N = 6: to layer 11). Two ingredients: the layer
  # is dual-unitary for any Z field (ZZ and the Y kick at π/4; ε only shifts the field), and
  # Z-basis product states are solvable initial states (from a random product state, or
  # with another ZZ angle, ⟨Z_j⟩ is not 0 at layer 1–2). So the KIM tests run t ≥ ⌈N/2⌉ + 1
  # layers with ε ≠ 0, or their ⟨Z_j⟩ checks would compare zeros.
  @testset "KIM shift-by-two invariance (dense, Néel), $name" for (name, circuit) in kim_sources
    # §2: with uniform ε the layer is invariant under j ↦ j mod N + 1 and, for even N, the
    # Néel state under a shift by two, so ⟨Z_j⟩ = ⟨Z_{j+2 mod N}⟩ at every layer.
    t = 6
    for N in (4, 6), ε in (0.1, 0.3)
      prof = kim_dense_profile(circuit(ε, t, N)..., N, kim_spec_L(N))
      @test size(prof) == (t + 1, N)
      shift2 = [mod1(j + 2, N) for j in 1:N]
      @test maximum(abs.(prof .- prof[:, shift2])) <= EXACT_ATOL
      # Not vacuous: the profile is not 0 after t = 0, and odd and even qubits differ at
      # some layer t ≥ 1 (so a profile uniform over the qubits, e.g. from |0…0⟩, is not
      # what passes here)
      @test maximum(abs, prof[2:end, :]) > 0.05
      @test maximum(abs.(prof[2:end, 1] .- prof[2:end, 2])) > 0.05
    end
  end

  @testset "KIM shift-by-two invariance needs the wrap bond (negative check)" begin
    # The open chain (the spec layer without (1, N)) is not shift invariant, so the check
    # above sees a dropped wrap bond.
    t = 6
    for N in (4, 6), ε in (0.1, 0.3)
      open_ = kim_spec_circuit(ε, t, N; wrap = false)
      prof = kim_dense_profile(open_..., N, kim_spec_L(N) - 1)
      shift2 = [mod1(j + 2, N) for j in 1:N]
      @test maximum(abs.(prof .- prof[:, shift2])) > 0.01
      # … and the wrap bond changes the profile
      @test maximum(abs.(prof .- kim_dense_profile(ε, t, N))) > 0.01
    end
  end

  @testset "KIM hand pins (dense), $name" for (name, circuit) in kim_sources
    # After layer 1 from the Néel state: ZZ and Z are diagonal, so the state is still the
    # basis state up to a phase; R_Y(π/2) = exp(−iπ/4·Y) then takes each Bloch vector
    # (0, 0, ±1) to (±1, 0, 0). So ⟨Z_j⟩ = 0 for every j, any N and any ε, and ⟨X_j⟩ = +1 on
    # odd j (|0⟩) and −1 on even j (|1⟩); |0…0⟩ or |1010…⟩ would give other signs.
    for N in 1:6, ε in (0.0, 0.1, -0.3, [0.05j - 0.2 for j in 1:N])
      strings, phases = circuit(ε, 1, N)
      ψ = DR.run_circuit(DR.basis_state(N, kim_neel(N)), strings, phases)
      b = [DR.bloch(ψ, j) for j in 1:N]
      @test all(abs(b[j][3]) <= EXACT_ATOL for j in 1:N)
      @test all(isapprox(b[j][1], isodd(j) ? 1.0 : -1.0; atol = EXACT_ATOL) for j in 1:N)
    end

    # N = 1 (no bond; Néel = |0⟩): a layer is R_Z(π/2 + ε) then R_Y(π/2) on the Bloch
    # vector. From (0, 0, 1): R_Z leaves it, R_Y gives (1, 0, 0), ⟨Z⟩ = 0. Layer 2: R_Z by
    # α = π/2 + ε gives (cos α, sin α, 0) = (−sin ε, cos ε, 0), and R_Y(π/2), z ↦ z cos − x sin,
    # gives z = sin ε. From |1⟩ every vector flips: −1, 0, −sin ε. These pin the reference
    # and N = 1, not the signs of the angles: with the Z phase flipped, layer 2 still gives
    # sin ε, and at N = 6 the ⟨Z_j⟩ profile is the same under a flipped ZZ, Z or Y sign, an X
    # kick or Z after Y. The structure test (symbols and phases at each position) and the
    # ⟨X_j⟩ pin above catch those.
    for ε in (0.0, 0.1, -0.3)
      strings, phases = circuit(ε, 2, 1)
      @test length(strings) == 4
      @test maxabsdiff(kim_dense_profile(strings, phases, 1, 2)[:, 1], [1, 0, sin(ε)]) <= EXACT_ATOL
      @test maxabsdiff(kim_dense_profile(strings, phases, 1, 2; ones = [1])[:, 1],
                       [-1, 0, -sin(ε)]) <= EXACT_ATOL
    end
  end
end
