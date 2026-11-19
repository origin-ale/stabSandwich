# Deprecated calls fail the suite (plan item 1). The suite's own process only warns on a
# deprecated call, and a warning is easy to miss, so this runs a typical disentangling path in
# a child process with `--depwarn=error`: any deprecated call on that path, in CampsPP or in a
# dependency it calls, throws there, the child exits non-zero and `run_child` throws.
#
# The path: `disentangler` on a string that takes every branch of `reducetoX`, `swapqubits`
# and `build_D` (all through `_apply_on_sites!`, where QuantumClifford 0.11.4 deprecated the
# old `apply!` argument order), then CAMPS-PP and CAMPS with the layer-end disentangler on a
# small XXZ circuit with general angles. Values are pinned elsewhere (TestDisentangler.jl,
# TestTimeEvolution.jl); here only that the path runs to the end.
#
# Slow tier only (STABSANDWICH_SLOW set): the child loads CampsPP from scratch, about 35 s.

if get(ENV, "STABSANDWICH_SLOW", "") != ""
@testset "no deprecated calls on a disentangling path (--depwarn=error)" begin
  N, t, ones = 6, 3, [2, 5]
  r = run_child("""
    using CampsPP, CampsPP.DisentangleCAMPS, QuantumClifford
    import MagicTensors as mt

    # Assumption pins: the flag reaches the child (2 = error) and turns a deprecation into an
    # error, here the old argument order apply!(stabilizer, op, sites) of QuantumClifford.
    report("depwarn option", Base.JLOptions().depwarn)
    report("old apply! order", errorinfo(() ->
      QuantumClifford.apply!(Stabilizer(one(CliffordOperator, 2).tab), tCNOT, [1, 2])))

    # k = 1, Q = I Z Y X Y: the first free X/Y is the Y on qubit 3, so a phase gate (reducetoX),
    # a swap 3 ↔ 2 (swapqubits), then controlled X, Y and Z from qubit 2 (build_D).
    D, sign = disentangler(1, P"IZYXY")
    report("disentangler", (nqubits(D), sign))

    gates, phases = xxz_circuit(0.3, 0.7, $t, $N)
    layer_ends = layerends($N, $t, xxz_circuit)
    ψ0 = computationalcamps(mt.QubitCAMPS, $N, $ones)
    O = transferredmagnetization($N, $ones)
    report("M", length(gates))
    evs, s, _, _, _, _, s_switch, _ = campspp_circuit_dynamics(
      ψ0, 4, 0.0, typemax(Int), gates, phases, O; layer_ends, track = true)
    report("campspp", (s, s_switch))
    _, i, _ = campssrc_circuit_dynamics(ψ0, gates, phases, 1e-11, O;
      criterion = :entangle, strategy = :brickwork, layer_ends)
    report("campssrc", i)
    """; flags = `--depwarn=error`)

  @test r["depwarn option"] == 2
  isdep, msg = r["old apply! order"]
  @test isdep && occursin("deprecated", msg)
  @test r["disentangler"] == (5, 1)   # D acts on the 5 qubits of Q; Q has sign +
  M = r["M"]
  s, s_switch = r["campspp"]
  @test s == M
  @test 0 < s_switch < M          # switched: both the CAMPS and the PP stage ran
  @test r["campssrc"] == M
end
end
