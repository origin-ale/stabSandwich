# Regression test that the PauliPropagation 0.7.3 hash bug is gone on PP ≥ 0.9 (the
# acceptance test of item 3 of the post-port plan), plus the tests of the Σ|c|² guard
# `check_pp_norm` (TimeEvolution.jl), which stays after the upgrade.
#
# The bug (PP 0.7.3, Julia ≥ 1.13): PP created the term types Julia lacks (UInt24, UInt40,
# UInt48, UInt56, …) with `@eval` the first time `getinttype` needed one. That added a `hash`
# method visible only from the next world age, so a `PauliSum` filled in the top-level
# statement that created its type could not find its own keys from the next statement on:
# terms were duplicated and Σ|c|² grew. PP ≥ 0.9 uses only built-in integers and
# `NTupleInteger` (CONVENTIONS.md §3), so nothing is created at runtime.
#
# Written before the upgrade: the bug checks fail on 0.7.3 by design (that red run is the
# proof that the trigger can see the bug) and must pass on ≥ 0.9.
#
# The bug checks run in fresh Julia processes that load only PauliPropagation (`run_child`,
# helpers.jl), so no workaround in CampsPP can mask the bug, and each N is the first use of
# its key type in that process. The children only report values; every check is made here.
# Expected values: hand values, the dense reference, the §3 key-type rule.

# -- The N that 0.7.3 broke (kept only to make sure the trigger covers them) ---------------

"N whose 0.7.3 key type (8⌈2N/8⌉ bits) was created at runtime with a BitIntegers `hash`."
const BROKEN_073_N = vcat(9:12, 17:28)
"The first N of each such 0.7.3 width (24, 40, 48, 56 bits): the statements that broke."
const FIRST_073_N = [9, 17, 21, 25]
"Both sides of every type boundary of the ≥ 0.9 rule, plus N = 46 and limb boundaries."
const RULE_BOUNDARY_N = [4, 5, 8, 9, 16, 17, 32, 33, 46, 64, 65, 96, 97, 128, 129]
"Trigger order: the 0.7.3 first-of-width N first, then 1:32, 33, 46 and the rule boundaries."
const TRIGGER_N = unique(vcat(FIRST_073_N, 1:32, 33, 46, RULE_BOUNDARY_N))

# -- Child code (`run_child` and its `report`/`errorinfo`/`datarows`: helpers.jl) ----------

"""
Child code: two top-level statements at N qubits. The first builds the one-term sum Z₁ (on
0.7.3 it created the term type if PP lacked it); the second adds Z₁ again. A correct sum then
holds one term with coefficient 2. Reports the sorted coefficients as "trigger N".
"""
trigger_code(N) = """
  P_$N = pp.PauliSum(pp.PauliString($N, [:Z], [1]))
  pp.add!(P_$N, pp.PauliString($N, [:Z], [1])); report("trigger $N", sort(collect(values(P_$N.terms))))
  """

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

@testset "PP hash bug gone (PauliPropagation ≥ 0.9) and the Σ|c|² guard" begin

  @testset "preconditions" begin
    # The bug showed only on Julia ≥ 1.13; on older Julia the trigger is clean even on 0.7.3,
    # so these tests would pass without proving anything.
    @test VERSION >= v"1.13"
    @test issubset(FIRST_073_N, BROKEN_073_N)
    @test issubset(BROKEN_073_N, TRIGGER_N)
    @test TRIGGER_N[1:length(FIRST_073_N)] == FIRST_073_N
    @test all(N -> N in TRIGGER_N, RULE_BOUNDARY_N)
  end

  @testset "bare PauliPropagation: the trigger keeps one term" begin
    r = run_child("import PauliPropagation as pp\n" * join(trigger_code.(TRIGGER_N)))
    @testset "N = $N" for N in TRIGGER_N
      @test r["trigger $N"] == [2.0]
    end
  end

  @testset "bare PauliPropagation: getinttype creates no types" begin
    # The root cause: on 0.7.3 `getinttype` evaluated new types into PP, which advances the
    # world counter. Up to 200 qubits (NTupleInteger{7}), nothing may be defined.
    r = run_child("""
      import PauliPropagation as pp
      let w = Base.get_world_counter()
        foreach(pp.getinttype, 1:200)
        report("world unchanged", Base.get_world_counter() == w)
      end
      """)
    @test r["world unchanged"]
  end

  @testset "bare PauliPropagation: the bug-report reproducer" begin
    # N = 10 rotation circuit, seed 101 (on 0.7.3 it gave 1.822 with Σc² = 4.34). The first
    # statement builds the circuit and Z₁, the first use of the key type in the child; the
    # second propagates. Gates are passed as literals, since the child does not load CampsPP.
    N = 10
    g, ph = rotation_circuit(MersenneTwister(101), 10, N)
    gates_code = "[" * join(("pp.PauliRotation($(repr(r.symbols)), $(repr(r.qinds)))" for r in g), ", ") * "]"
    r = run_child("""
      import PauliPropagation as pp
      g = $gates_code; ph = $(repr(ph)); Z1 = pp.PauliSum(pp.PauliString($N, [:Z], [1]))
      S = pp.propagate(g, Z1, -2 .* ph; min_abs_coeff = 0.0); report("reproducer", (pp.overlapwithzero(S), sum(abs2, values(S.terms))))
      """)
    Z1 = pp.PauliSum(pp.PauliString(N, [:Z], [1]))
    ref = dense_dynamics(N, Int[], g, ph, Z1, [length(g)])[end]
    ev, sumsq = r["reproducer"]
    @test abs(ev - ref) <= EXACT_ATOL
    @test abs(sumsq - 1) <= EXACT_ATOL
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

  # Not kept from the 0.7.3 file: "pauliprop_circuit_dynamics throws on a corrupt observable",
  # which raised Σ|c|² with the hash bug itself. With the bug gone nothing in `propagate` can
  # raise it: Pauli rotations and Cliffords preserve it, PP rejects noise strengths outside
  # [0, 1] and has no Heisenberg-picture ImaginaryPauliRotation (checked on 0.7.3 and 0.9.2),
  # and a custom gate would need PP internals that changed between the versions. A sum with
  # an inflated Σ|c|² built by hand does not trip the guard either, since the guard compares
  # against the input's own Σ|c|². So only the in-process `check_pp_norm` tests remain.
end
