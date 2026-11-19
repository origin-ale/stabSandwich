# Thread safety (plan item 2). TransferredMag_varparams.jl runs samples in `Threads.@threads`;
# at -t 2 it used to fail with "rotation needs a Hermitian Pauli string". The cause: in
# QuantumClifford 0.11, `apply!(stab, op[, sites])` and `l * r` use `op.buffer` (`l.buffer`) as
# scratch space, and CampsPP passed shared constants there (QC's `tSWAP`, `tPhase`, `tCNOT`
# and the `_OLD_R` cores of `old_entangle_gates`), so concurrent tasks overwrote each other's
# rows and phases.
#
# Spec: a run gives the same result whatever else runs at the same time. Reference: the same
# seeds run serially (plan, "How tests are written"). Three levels:
#   1. fast tier, in this process: no shared `CliffordOperator` constant has its buffer
#      written on any path a threaded script reaches. This is the precondition for the race,
#      so it catches a regression deterministically, without threads. With the QC pins it
#      relies on (`apply!` and `*` write the buffer, and a `copy` has its own).
#   2. slow tier, child process at -t 4: stress tests of `_apply_on_sites!`, `disentangler`
#      and `old_entangle_gates`, and the dynamics paths, each compared bitwise with a serial
#      run of the same seeds in the same process (so BLAS threading is the same for both).
#      Bitwise equality was checked to hold before it was required (2026-11-19).
#   3. slow tier, root test/TestThreadedScript.jl: the script itself at -t 1, 2 and 4.
# Sensitivity (2026-11-19, pre-fix code patched in memory in scratch processes, src untouched;
# mutant A = `_apply_on_sites!` without the copy, mutant B = `old_entangle_gates` without it):
#   - level 1: 5 of the 13 "no shared constant is written" checks fail under A and 5 under B,
#     deterministically; the pins and the scan pass under both;
#   - level 2 (3 child runs per mutant, all failed): under A, `_apply_on_sites!` and
#     `disentangler` mismatch and `camps` errors with "rotation needs a Hermitian Pauli string"
#     in every run, `campspp` in 2 of 3; under B, `old_entangle_gates` and
#     `BondSearchDisentangler` mismatch in every run, `campssrc entangle/brickwork` errors.

const DC = CampsPP.DisentangleCAMPS

# -- Shared Clifford constants ------------------------------------------------------------

"""
Every `CliffordOperator` held by a constant global of the modules `mods`, directly or in a
tuple, as `"Module.name" => op`; each object once (`_OLD_R[1]` is `tSWAP`).
"""
function shared_cliffords(mods...)
  found = Pair{String, CliffordOperator}[]
  for m in mods, n in names(m; all = true)
    (isdefined(m, n) && isconst(m, n)) || continue
    v = getglobal(m, n)
    if v isa CliffordOperator
      push!(found, "$m.$n" => v)
    elseif v isa Tuple
      for (i, x) in enumerate(v)
        x isa CliffordOperator && push!(found, "$m.$n[$i]" => x)
      end
    end
  end
  return unique(p -> objectid(p.second), found)
end

"""
Run `f()` with the scratch buffer of every operator in `ops` (as from `shared_cliffords`) set
to a sentinel, and return the names of those whose buffer `f` changed. The buffers are
restored afterwards (QC clears a buffer before each use, so their contents don't matter).
"""
function written_buffers(f, ops)
  saved = [(copy(op.buffer.xz), op.buffer.phase[]) for (_, op) in ops]
  sentinel = [rand(MersenneTwister(i), UInt64, length(op.buffer.xz)) for (i, (_, op)) in enumerate(ops)]
  try
    for ((_, op), s) in zip(ops, sentinel)
      op.buffer.xz .= s
      op.buffer.phase[] = 0x3
    end
    f()
    return [name for ((name, op), s) in zip(ops, sentinel)
            if op.buffer.xz != s || op.buffer.phase[] != 0x3]
  finally
    for ((_, op), (xz, ph)) in zip(ops, saved)
      op.buffer.xz .= xz
      op.buffer.phase[] = ph
    end
  end
end

"A doped 8-qubit XXZ sample as in TransferredMag_varparams.jl: (gates, phases, ψ, ones, O, layer ends)."
function threads_sample(seed; N = 8, t = 4, p = 0.5)
  rng = MersenneTwister(seed)
  gates, phases = xxz_circuit(π/4, π/4, t, N)
  phases = xy_magic(rng, phases, p; magicphase = 3π/16)
  ψ, ones = domainwallstate(rng, N, 0.6)
  return gates, phases, ψ, ones, transferredmagnetization(N, ones), layerends(N, t, xxz_circuit)
end

@testset "thread safety (plan item 2)" begin
  ops = shared_cliffords(QuantumClifford, DC, CampsPP, mt)

  @testset "shared constants found" begin
    # The scan is not vacuous: it holds the constants the pre-fix code wrote to.
    for x in (tSWAP, tPhase, tCNOT, DC._OLD_R...)
      @test any(op === x for (_, op) in ops)
    end
  end

  @testset "QuantumClifford: which buffer is scratch (assumption pins)" begin
    # These pins are also the negative checks of `written_buffers`: a write is detected.
    C = one(CliffordOperator, 3)
    name(x) = only(n for (n, op) in ops if op === x)
    @test written_buffers(() -> apply!(Stabilizer(C.tab), [1, 3], tSWAP), ops) == [name(tSWAP)]
    @test written_buffers(() -> apply!(Stabilizer(C.tab), [2], tPhase), ops) == [name(tPhase)]
    @test written_buffers(() -> tCNOT * tSWAP, ops) == [name(tCNOT)]          # l * r: l's buffer
    @test written_buffers(() -> tId1 ⊗ tPhase, ops) == String[]
    # The fix relies on this: a copy equals the original and has its own buffer, so applying
    # it, or multiplying by it, leaves the original untouched.
    for (_, op) in ops
      c = copy(op)
      @test c == op
      @test c.buffer !== op.buffer && c.buffer.xz !== op.buffer.xz && c.buffer.phase !== op.buffer.phase
    end
    @test isempty(written_buffers(() -> apply!(Stabilizer(C.tab), [1, 3], copy(tSWAP)), ops))
    @test isempty(written_buffers(() -> copy(tCNOT) * tSWAP, ops))
  end

  @testset "no shared constant is written" begin
    # Fails deterministically on the pre-fix code (checked 2026-11-19): `_apply_on_sites!`
    # wrote tSWAP, tPhase and tCNOT; `old_entangle_gates` wrote the `_OLD_R` cores.
    # disentangler: Y on the chosen qubit, a swap, and controlled X, Y, Z (as TestDepwarn.jl);
    # then X already on qubit k+1.
    @test isempty(written_buffers(() -> disentangler(1, P"IZYXY"), ops))
    @test isempty(written_buffers(() -> disentangler(0, P"-XZY"), ops))
    rng = MersenneTwister(7)
    C = random_clifford(rng, 5)
    @test isempty(written_buffers(() -> disentangler(2, C, P"YXZYX"), ops))
    @test isempty(written_buffers(() -> for (op, s) in ((tSWAP, [1, 4]), (tPhase, [2]), (tCNOT, [5, 3]))
      DC._apply_on_sites!(C, op, s)
    end, ops))
    @test isempty(written_buffers(DC.old_entangle_gates, ops))
    @test isempty(written_buffers(() -> DC.BondSearchDisentangler(:chi3), ops))
    @test isempty(written_buffers(() -> sweep_disentangler(:entangle, :brickwork, 8), ops))

    # Every dynamics path a threaded script can reach, on one doped sample.
    gates, phases, ψ, ones, O, lends = threads_sample(3)
    M = length(gates)
    local r
    @test isempty(written_buffers(() -> (r = campspp_circuit_dynamics(ψ, 4, 1e-6, 1_000_000,
      gates, phases, O; layer_ends = lends, track = true)), ops))
    @test 0 < r[7] < M && r[2] == M          # switched: CAMPS, the seam and PP all ran
    @test isempty(written_buffers(() -> camps_circuit_dynamics(ψ, gates, phases, typemax(Int), O;
      layer_ends = lends), ops))
    @test isempty(written_buffers(() -> pauliprop_circuit_dynamics(ones, gates, phases, 1e-10,
      1_000_000, O; layer_ends = lends), ops))
    for (criterion, strategy) in ((:chi3, :snake), (:entangle, :brickwork))
      @test isempty(written_buffers(() -> campssrc_circuit_dynamics(ψ, gates, phases, 1e-10, O;
        criterion, strategy, layer_ends = lends), ops))
    end
  end
end

# -- Concurrent runs (child process at -t 4) ----------------------------------------------

"""
Child code: each case run for seeds 1:4nthreads serially, then several times with one task per
seed, all spawned at once; every concurrent result compared bitwise (`isequal`) with the serial
one and checked with the case's `valid`. Reports per case a NamedTuple (see `stress`).
"""
const THREADS_CHILD = raw"""
using CampsPP, CampsPP.DisentangleCAMPS, QuantumClifford, Random
using LinearAlgebra: BLAS
const DC = CampsPP.DisentangleCAMPS
BLAS.set_num_threads(Threads.nthreads())      # as TransferredMag_varparams.jl (setup_threads)
report("threads", (Threads.nthreads(), BLAS.get_num_threads()))

"Every image in the tableau has a real sign (phase +1 or −1), i.e. is Hermitian."
hermitian(C::CliffordOperator) = all(iseven, QuantumClifford.phases(C))

function stress(f, name; valid = _ -> true, trials, ntask = 4 * Threads.nthreads())
  ref = [f(s) for s in 1:ntask]
  nbad = nerr = ninvalid = 0
  firsterr = ""
  tids = Set{Int}()
  lk = ReentrantLock()
  for _ in 1:trials
    tasks = [Threads.@spawn(begin
               lock(() -> push!(tids, Threads.threadid()), lk)
               f(s)
             end) for s in 1:ntask]
    for (s, task) in enumerate(tasks)
      try
        res = fetch(task)
        isequal(res, ref[s]) || (nbad += 1)
        valid(res) || (ninvalid += 1)
      catch e
        nerr += 1
        e isa TaskFailedException && (e = e.task.exception)
        isempty(firsterr) && (firsterr = first(sprint(showerror, e), 300))
      end
    end
  end
  report(name, (; nbad, nerr, ninvalid, nruns = trials * ntask, nthreads_used = length(tids),
                refvalid = all(valid, ref), firsterr))
end

# _apply_on_sites! with QC's shared gates, 200 random gates on a random 6-qubit frame
function apply_case(seed)
  rng = MersenneTwister(seed)
  N = 6
  C = random_clifford(rng, N)
  for _ in 1:200
    g = rand(rng, 1:3)
    g == 2 ? DC._apply_on_sites!(C, tPhase, [rand(rng, 1:N)]) :
             DC._apply_on_sites!(C, g == 1 ? tSWAP : tCNOT, randperm(rng, N)[1:2])
  end
  return C
end

# disentangler on random disentanglable (k, C, P); each result with CONVENTIONS.md §4
# property 1, D†QD = s·X_{k+1} (C†PC as cliffconj(inv(C), P), §4)
cliffconj(C, P) = (C * Stabilizer([P]))[1]
function disentangler_case(seed)
  rng = MersenneTwister(seed)
  N = 6
  out = []
  while length(out) < 40
    C = random_clifford(rng, N)
    k = rand(rng, 0:N-1)
    P = random_pauli(rng, N; nophase = false, realphase = true)
    Q = cliffconj(inv(C), P)
    any(xbit(Q)[k+1:end]) || continue
    D, s = disentangler(k, C, P)
    X = zero(PauliOperator, N); X[k+1] = (true, false)
    push!(out, (D, s, cliffconj(inv(D), Q) == s * X))
  end
  return out
end
disentangler_valid(out) = all(((D, s, prop),) -> hermitian(D) && s in (1, -1) && prop, out)

old_case(seed) = [DC.old_entangle_gates() for _ in 1:10]
old_valid(v) = all(gs -> length(gs) == 19 && all(hermitian, gs), v)
bond_case(seed) = DC.BondSearchDisentangler(:chi3).matrices

stress(apply_case, "_apply_on_sites!"; valid = hermitian, trials = 5)
stress(disentangler_case, "disentangler"; valid = disentangler_valid, trials = 5)
stress(old_case, "old_entangle_gates"; valid = old_valid, trials = 5)
stress(bond_case, "BondSearchDisentangler"; trials = 5)

# Dynamics paths, one doped 8-qubit sample per seed, as in TransferredMag_varparams.jl
function sample(seed; N = 8, t = 4, p = 0.5)
  rng = MersenneTwister(seed)
  gates, phases = xxz_circuit(π/4, π/4, t, N)
  phases = xy_magic(rng, phases, p; magicphase = 3π/16)
  ψ, ones = domainwallstate(rng, N, 0.6)
  return gates, phases, ψ, ones, transferredmagnetization(N, ones), layerends(N, t, xxz_circuit)
end
const M = 84                                  # gates in 4 layers at N = 8
finite(evs) = all(isfinite, evs)

function campspp_case(seed)
  gates, phases, ψ, _, O, lends = sample(seed)
  evs, s, bonddims, nterms, _, _, s_switch, k_switch = campspp_circuit_dynamics(
    ψ, 4, 1e-6, 1_000_000, gates, phases, O; layer_ends = lends, track = true)
  return (evs, s, bonddims, nterms, s_switch, k_switch)
end
campspp_valid(r) = finite(r[1]) && r[2] == M && 0 < r[5] < M     # a real switch
function camps_case(seed)
  gates, phases, ψ, _, O, lends = sample(seed)
  _, k, i, evs, bonddims = camps_circuit_dynamics(ψ, gates, phases, typemax(Int), O;
    layer_ends = lends, track = true)
  return (evs, k, i, bonddims)
end
function pp_case(seed)
  gates, phases, _, ones, O, lends = sample(seed)
  return pauliprop_circuit_dynamics(ones, gates, phases, 1e-10, 1_000_000, O;
    layer_ends = lends, track = true)
end
function campssrc_case(seed, criterion, strategy)
  gates, phases, ψ, _, O, lends = sample(seed)
  _, i, evs, bonddims = campssrc_circuit_dynamics(ψ, gates, phases, 1e-10, O;
    criterion, strategy, layer_ends = lends, track = true)
  return (evs, i, bonddims)
end

stress(campspp_case, "campspp"; valid = campspp_valid, trials = 2)
stress(camps_case, "camps"; valid = r -> finite(r[1]) && r[3] == M, trials = 2)
stress(pp_case, "pauliprop"; valid = r -> finite(r[2]) && r[1] == M, trials = 2)
stress(s -> campssrc_case(s, :chi3, :snake), "campssrc chi3/snake";
  valid = r -> finite(r[1]) && r[2] == M, trials = 2)
stress(s -> campssrc_case(s, :entangle, :brickwork), "campssrc entangle/brickwork";
  valid = r -> finite(r[1]) && r[2] == M, trials = 2)
"""

const THREADS_CASES = ["_apply_on_sites!", "disentangler", "old_entangle_gates",
  "BondSearchDisentangler", "campspp", "camps", "pauliprop", "campssrc chi3/snake",
  "campssrc entangle/brickwork"]

"The checks on the child's reports (also used by the sensitivity runs)."
function check_threads_child(r)
  @test r["threads"] == (4, 4)
  @testset "$name" for name in THREADS_CASES
    c = r[name]
    @test c.refvalid                       # the serial reference itself passes `valid`
    @test c.nthreads_used >= 2             # the tasks really ran concurrently
    @test c.nerr == 0
    c.nerr == 0 || @info "$name: first error" c.firsterr
    @test c.nbad == 0                      # bitwise equal to the serial run
    @test c.ninvalid == 0
  end
end

# Slow tier only: the child loads CampsPP from scratch (about 35 s) and compiles every path.
if get(ENV, "STABSANDWICH_SLOW", "") != ""
@testset "concurrent runs equal serial runs (child, -t 4)" begin
  check_threads_child(run_child(THREADS_CHILD; flags = `-t 4`))
end
end
