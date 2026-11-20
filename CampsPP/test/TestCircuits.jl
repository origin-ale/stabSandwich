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
end
