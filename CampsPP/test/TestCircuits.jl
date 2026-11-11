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
