# The PauliPropagation 0.7.3 hash bug and CampsPP's two fixes:
#   1. `CampsPP.__init__` (CampsPP.jl) creates every PP term type up to 64 bits (N ≤ 32)
#      when CampsPP loads;
#   2. `check_pp_norm` (TimeEvolution.jl), called after each `pp.propagate` in
#      `pauliprop_circuit_dynamics`, throws when Σ|c|² rises above the input observable's.
#
# The bug: PP creates the term types Julia lacks (UInt24, UInt40, UInt48, UInt56 below 64
# bits) with `@eval` the first time `getinttype` needs one. That adds a `hash` method visible
# only from the next world age, so a `PauliSum` filled in the top-level statement that created
# its type cannot find its own keys from the next statement on: terms are duplicated and
# Σ|c|² grows. Affected: N = 9–12 and 17–28. Diagnosis: port plan, "Pre-existing upstream bug
# found in the smoke runs".
#
# Whether a type exists depends on everything a process has run, and this process loaded
# CampsPP long ago, so the load-order checks run in fresh Julia processes on the test
# environment (`run_child`). The children only report values; every check is made here.
# Expected values: the CONVENTIONS.md §3 key-width rule, hand values, the dense reference.
#
# On upgrading to PauliPropagation ≥ 0.9 (no types created at runtime) the bare-PP pin and the
# type checks of "after `using CampsPP`" fail by design: delete them with the `__init__`
# workaround, and keep the Σ|c|² guard tests.

# -- Key widths ---------------------------------------------------------------------------

"Key width in bits of a PP term at N qubits (CONVENTIONS.md §3): 8⌈2N/8⌉."
keybits(N) = 8 * cld(2N, 8)

# Widths below 64 bits that Julia has no built-in type for: PP creates them at runtime, and
# BitIntegers gives them their own `hash`.
const RUNTIME_WIDTHS = [24, 40, 48, 56]
const AFFECTED_N = [N for N in 1:32 if keybits(N) in RUNTIME_WIDTHS]
"The first N of each runtime width, whose first use creates the type."
const FIRST_N = [minimum(N for N in AFFECTED_N if keybits(N) == w) for w in RUNTIME_WIDTHS]

# -- Child code (`run_child` and its `report`/`errorinfo`/`datarows`: helpers.jl) ----------

"""
Child code: two top-level statements at N qubits. The first builds the one-term sum Z₁, which
creates the term type if PP lacks it; the second adds Z₁ again. A correct sum then holds one
term with coefficient 2. Reports the sorted coefficients as "trigger N".
"""
trigger_code(N) = """
  P_$N = pp.PauliSum(pp.PauliString($N, [:Z], [1]))
  pp.add!(P_$N, pp.PauliString($N, [:Z], [1])); report("trigger $N", sort(collect(values(P_$N.terms))))
  """

"Child code: report whether PP defines each runtime-width type, as `name`."
defined_code(name) =
  "report(\"$name\", [isdefined(pp, Symbol(:UInt, w)) for w in $RUNTIME_WIDTHS])\n"

# -- Guard helpers ------------------------------------------------------------------------

const GUARD_STRINGS = [([:Z], [1]), ([:Z], [2]), ([:Z], [3]), ([:X], [1])]

"A 3-qubit pp sum whose terms are Z₁, Z₂, Z₃, X₁ (in that order) with coefficients `cs`."
function sum_with(cs; T = Float64)
  S = pp.PauliSum(T, 3)
  for (c, (syms, qs)) in zip(cs, GUARD_STRINGS)
    pp.add!(S, pp.PauliString(3, syms, qs, T(c)))
  end
  return S
end

"A one-term sum with Σ|c|² = `x` (up to rounding)."
sum_of_squares(x) = sum_with([sqrt(x)])

@testset "PP hash bug (PauliPropagation 0.7.3) and its fixes" begin

  @testset "runtime key widths" begin
    @test AFFECTED_N == vcat(9:12, 17:28)
    @test FIRST_N == [9, 17, 21, 25]
    # PP names its types by the §3 width, so the names the children look up are the right ones.
    @test all(nameof(pp.getinttype(N)) == Symbol(:UInt, keybits(N)) for N in [1:32; 46])
  end

  @testset "bare PauliPropagation: the trigger corrupts the sum (assumption pin)" begin
    # Without CampsPP, the first use of each runtime width creates the type in the statement
    # that fills the sum, and the next statement cannot find Z₁, so it adds a second entry
    # beside the stale one. Unaffected: N = 12, whose type N = 9 already made; N = 14
    # (built-in UInt32); N = 33 and 46 (UInt72, UInt96: no BitIntegers `hash` above 64 bits;
    # N = 46 is the production size). This shows the trigger can see the bug, so the
    # fixed-side test below is not vacuous.
    others = [12, 14, 33, 46]
    r = run_child("import PauliPropagation as pp\n" * defined_code("defined before") *
                  join(trigger_code.(vcat(FIRST_N, others))))
    @test r["defined before"] == falses(length(RUNTIME_WIDTHS))
    @testset "N = $N, creates UInt$(keybits(N))" for N in FIRST_N
      @test r["trigger $N"] == [1.0, 1.0]
    end
    @testset "N = $N, unaffected" for N in others
      @test r["trigger $N"] == [2.0]
    end
  end

  @testset "after `using CampsPP`" begin
    # Fix 1: loading CampsPP creates the runtime-width types, so nothing at N ≤ 32 creates a
    # type later. The reproducer from the bug report (its value was 1.822 with Σc² = 4.34)
    # runs before the triggers, so it is the first use of UInt24 if `__init__` did nothing.
    r = run_child("import PauliPropagation as pp\n" * defined_code("defined before") *
                  "using CampsPP\n" * defined_code("defined after") * """
      let w = Base.get_world_counter()
        foreach(pp.getinttype, 1:32)
        report("world unchanged", Base.get_world_counter() == w)
      end
      using Random
      g, ph = rotation_circuit(MersenneTwister(101), 10, 10)
      S = pp.propagate(g, pp.PauliSum(pp.PauliString(10, [:Z], [1])), -2 .* ph; min_abs_coeff = 0.0)
      report("reproducer", (pp.overlapwithzero(S), sum(abs2, values(S.terms))))
      """ * join(trigger_code.(AFFECTED_N)))
    @test r["defined before"] == falses(length(RUNTIME_WIDTHS))  # so CampsPP made them
    @test r["defined after"] == trues(length(RUNTIME_WIDTHS))
    @test r["world unchanged"]           # `getinttype` defines nothing for N ≤ 32 any more

    g, ph = rotation_circuit(MersenneTwister(101), 10, 10)
    Z1 = pp.PauliSum(pp.PauliString(10, [:Z], [1]))
    ref = dense_dynamics(10, Int[], g, ph, Z1, [length(g)])[end]
    ev, sumsq = r["reproducer"]
    @test abs(ev - ref) <= EXACT_ATOL
    @test abs(sumsq - 1) <= EXACT_ATOL

    @testset "trigger, N = $N" for N in AFFECTED_N
      @test r["trigger $N"] == [2.0]
    end
  end

  @testset "sumsquares and check_pp_norm" begin
    sumsq = CampsPP.sumsquares
    check = CampsPP.check_pp_norm

    @test sumsq(sum_with([0.5, -0.5, 0.5, -0.5])) === 1.0
    @test sumsq(sum_with([0.6im, 0.8]; T = ComplexF64)) isa Float64
    @test abs(sumsq(sum_with([0.6im, 0.8]; T = ComplexF64)) - 1) <= 1e-15
    @test sumsq(pp.PauliSum(3)) === 0.0

    # Σ|c|² of the output against that of the input observable, here exactly 1.
    @test check(1.0, sum_with([0.5, 0.5, 0.5, 0.5])) === nothing   # equal
    @test check(1.0, sum_with([0.5, 0.5])) === nothing             # truncated
    @test check(1.0, pp.PauliSum(3)) === nothing                    # everything truncated
    @test check(1.0, sum_with([0.6im, 0.8]; T = ComplexF64)) === nothing
    @test check(1.0, sum_of_squares(1 + 5e-9)) === nothing          # rounding headroom
    @test_throws ErrorException check(1.0, sum_of_squares(1 + 2e-8))
    @test_throws ErrorException check(1.0, sum_with([0.6im, 0.9]; T = ComplexF64))
    @test_throws "corrupt" check(1.0, sum_with([1.0, 1.0]))

    # The tolerance is relative: each of these excesses would be judged the other way by an
    # absolute 1e-8.
    @test check(100.0, sum_of_squares(100 * (1 + 5e-9))) === nothing   # excess 5e-7
    @test_throws ErrorException check(100.0, sum_of_squares(100 * (1 + 2e-8)))
    @test_throws ErrorException check(0.01, sum_of_squares(0.01 * (1 + 5e-8)))  # excess 5e-10

    @test check(1.0, sum_of_squares(1.5); rtol = 0.6) === nothing
    @test_throws ErrorException check(1.0, sum_of_squares(1.5); rtol = 0.4)
  end

  @testset "untruncated propagation keeps Σ|c|² far inside the tolerance" begin
    # The premise of the guard: rotations preserve Σ|c|², so only rounding moves it, by much
    # less than the relative 1e-8 that `check_pp_norm` allows.
    N, t = 6, 4
    circuits = [xxz_circuit(0.3, 0.7, t, N), xxz_circuit(π/4, 3π/16, t, N),
                rotation_circuit(Xoshiro(5), 40, N)]
    for (gates, phases) in circuits, O in (transferredmagnetization(N, [1, 2, 3]),
                                           pp.PauliSum(pp.PauliString(N, [:Z], [1])))
      norm0 = sum(abs2, values(O.terms))
      drift = maximum(eachindex(gates)) do l
        S = pp.propagate(gates[1:l], O, -2 .* phases[1:l]; min_abs_coeff = 0.0)
        abs(sum(abs2, values(S.terms)) / norm0 - 1)
      end
      @test drift <= 1e-12
    end
  end

  @testset "pauliprop_circuit_dynamics throws on a corrupt observable" begin
    # Fix 2. The child builds Z₁ at N = 10 with bare PP, in the statement that creates UInt24,
    # and loads CampsPP after it, so the observable is corrupt (checked) and the first
    # propagation raises Σ|c|². Every PP entry point must throw at the first record point,
    # before writing its value.
    N, ones, t, χ_campspp = 10, [2, 5], 2, 2
    le = spec_layer_ends(N, t)
    r = run_child("""
      import PauliPropagation as pp
      obs = pp.PauliSum(pp.PauliString($N, [:Z], [1]))
      report("obs corrupt", let C = deepcopy(obs)
        pp.add!(C, pp.PauliString($N, [:Z], [1])); sort(collect(values(C.terms)))
      end)
      using CampsPP
      import MagicTensors as mt
      gates, phases = xxz_circuit(0.3, 0.7, $t, $N)
      le = $le
      ψ = computationalcamps(mt.QubitCAMPS, $N, $ones)
      mktempdir() do dir
        path = joinpath(dir, "onebit.txt")
        initialize_output(path, "Z1", Dict())
        report("onebit", errorinfo(() -> pauliprop_circuit_dynamics(
          $ones, gates, phases, 0.0, typemax(Int), obs, path; layer_ends = le)))
        report("onebit rows", datarows(path))
        path = joinpath(dir, "camps.txt")
        initialize_output(path, "Z1", Dict())
        report("camps", errorinfo(() -> pauliprop_circuit_dynamics(
          ψ, 0, gates, phases, 0.0, typemax(Int), obs, path; layer_ends = le)))
        report("camps rows", datarows(path))
      end
      report("campspp", errorinfo(() -> campspp_circuit_dynamics(
        ψ, $χ_campspp, 0.0, typemax(Int), gates, phases, obs; layer_ends = le)))
      """)
    @test r["obs corrupt"] == [1.0, 1.0]   # the bug did hit the observable
    @testset "$entry" for entry in ("onebit", "camps", "campspp")
      threw, msg = r[entry]
      @test threw
      @test occursin("corrupt", msg)
    end
    @test r["onebit rows"] == [0]          # the t = 0 value only
    @test r["camps rows"] == Int[]         # this method records no t = 0 value

    @testset "control: the same runs on a correct observable" begin
      gates, phases = xxz_circuit(0.3, 0.7, t, N)
      O = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
      ref = dense_dynamics(N, ones, gates, phases, O, le)
      _, evs = pauliprop_circuit_dynamics(ones, gates, phases, 0.0, typemax(Int), O; layer_ends = le)
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      ψ = computationalcamps(mt.QubitCAMPS, N, ones)
      _, evs = pauliprop_circuit_dynamics(ψ, 0, gates, phases, 0.0, typemax(Int), O; layer_ends = le)
      @test maxabsdiff(evs, ref[2:end]) <= EXACT_ATOL
      evs, _, _, _, _, _, s_switch = campspp_circuit_dynamics(
        ψ, χ_campspp, 0.0, typemax(Int), gates, phases, O; layer_ends = le, track = true)
      @test maxabsdiff(evs, ref) <= EXACT_ATOL
      @test s_switch < le[1]                 # PP evaluates every record point
    end
  end
end
