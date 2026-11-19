# Convention spec (T2)

**Status: approved by the user on 2026-11-11** (five decisions, listed at the end). Revised the same day after two independent ambiguity reviews; the revision changed wording and added definitions, not conventions (see the change log). **Updated 2026-11-13 (user approval) for the port:** the new frame (§4) and the new converter (§3) are no longer provisional, and the Clifford-angle tolerance of the port is fixed (§2). **Updated 2026-11-16 (user approval) after the T9 review:** invalid-input handling, keyword forwarding and the layer-end disentangler (§2); MPS renormalisation and the scope of the truncated-regime check (§6). **Updated 2026-11-17 (user approval) after Phase 6 deleted the old code:** the old-code items are kept, marked historical, since the saved old data and the T11 comparisons refer to them; QuantumClifford has no `C * PauliOperator` (§4); one entropy function family (§5). **Updated 2026-11-19 (user approval) for PauliPropagation ≥ 0.9:** the key types (§3); and for post-port item 4: the pp seam (§5).

Every expected value in the tests comes from this page, the dense reference or the hand-derived golden set, never from the implementation.

**Terms used on this page**
- *Old code*: the pre-port packages (`DisentangleCAMPS`, now the submodule `CampsPP.DisentangleCAMPS`, and `CampsPP`, on the vendored `CliffordMPS`), with the state type `CliffordMPS.CAMPS`. *New code*: the port to MagicTensors.jl (pinned at commit `72eac4b4310042ee714d784e4762121b4f965b68`), with the state type `MagicTensors.QubitCAMPS`. Both stacks lived side by side until Phase 6, which deleted the old code (2026-11-17); it stays recoverable from the `pre-port` tag, and its items below are historical.
- *Dense reference*: `CampsPP/test/DisentangleCAMPS/DenseReference.jl`, an exact state-vector simulator written without reading the implementation.
- *Golden set / assumption pins*: hand-derived checks in `CampsPP/test/DisentangleCAMPS/TestGoldenSet.jl` and `TestEvolve.jl`.
- *PP*: PauliPropagation.jl (≥ 0.9; 0.9.2 when this was written). *QC*: QuantumClifford.jl (0.11.4).
- Symbols: ϕ is always a gate angle in the exp(iϕP) convention; θ is only ever a PP angle; α is the angle of a magic state; |φ⟩ is an arbitrary state.

## 1. Qubits, basis, states

- Qubits are numbered 1…N. Qubit j is MPS site j, position j of a QC string (`P"XYZ"` = X₁Y₂Z₃), and character j of a saved bitstring.
- |0⟩ is the +1 eigenstate of Z, |1⟩ the −1 eigenstate. Paulis are the standard matrices, with Y = iXZ = [[0, −i], [i, 0]].
- `onebitinds` (also called `init_ones`) is the list of qubits in |1⟩. `computationalcamps(N, onebitinds)` is that basis state.
- **Domain-wall ensemble.** `domainwallstate(rng, N, μ)` returns a random basis state and its `onebitinds`. With p = e^μ / (e^μ + e^−μ), each qubit of the left half (1…⌊N/2⌋) is |1⟩ independently with probability 1 − p, and each qubit of the right half with probability p. For μ > 0 the excitations start mostly on the right; μ = 0 gives probability 1/2 everywhere.
- Dense reference: the amplitude of |b₁…b_N⟩ is stored at index 1 + Σⱼ bⱼ 2^(j−1), so qubit 1 is the least significant bit and A₁⊗…⊗A_N is `kron(A_N, …, A_1)`. Comparisons with the code go through Pauli expectation values only, so no memory layout is shared.

## 2. Gates and circuits

- **Rotation.** A gate is a pair (P, ϕ), with P a Hermitian Pauli string (sign ±1), and means U = exp(iϕP) = cos ϕ · 1 + i sin ϕ · P for every real ϕ, sign included.
  - Special values: exp(±iπ/4·P) = (1 ± iP)/√2, exp(±iπ/2·P) = ±iP, exp(i·3π/4·P) = (−1 + iP)/√2, exp(−i·3π/4·P) = −exp(iπ/4·P), and ϕ = 0, π give ±1.
  - exp(iπ/4·P) maps A ↦ iPA for every Pauli A that anticommutes with P. In particular exp(iπ/4·Z) is S† up to a phase, and the phase gate S = diag(1, i) is exp(−iπ/4·Z) up to a phase.
  - **Signed strings.** exp(iϕ(−P)) = exp(−iϕP).
  - **Invalid input.** The new code rejects invalid arguments with `ArgumentError`, the one exception type used for them: a gate string with phase ±i (not Hermitian, so exp(iϕP) is not unitary), in both `apply!` methods and in `rotation_gate`; and a non-finite angle (NaN, ±Inf) in both `apply!` methods.
- **Known bugs in the old `apply!`** (both methods; historical). These were the only `@test_broken` cases, and only on the old type; the new code passes them, and since Phase 6 they are plain tests:
  1. ϕ = −π/4 and ϕ = +3π/4 are wrong in every frame, including the identity frame;
  2. ϕ = ±π/2 is wrong when the frame does not map P to ±P, i.e. when C†PC ≠ ±P (the old code multiplies the frame by P on the wrong side). This can only happen once the frame holds a non-Pauli Clifford, i.e. after a disentangler or a π/4 gate; with a Pauli frame, or whenever C†PC = ±P, the branch is correct.
- **Clifford angles and tolerance.** The physical state must equal exp(iϕP)|ψ⟩ for every real ϕ. Which angles an implementation recognises as Clifford affects only frame-dependent quantities (k, the frame, bond dimensions).
  - The old code recognises |ϕ| ∈ {0, π/4, π/2, 3π/4, π}; other multiples of π/4 (5π/4, 2π, …) take the generic branch. It uses `atol = 1e-6` at 0 and the default `≈` relative tolerance elsewhere.
  - The new code makes this uniform (user, 2026-11-11; defaults decided 2026-11-12): ϕ is taken as Clifford when `isapprox(ϕ, mπ/4; atol, rtol)` holds for the nearest integer m, for every m (any multiple of π/4, of any size, including 0 and ±π). Both `apply!` methods (the four-argument one and the truncating `apply!(ψ, P, ϕ, thl)`) take `atol` and `rtol` as keyword arguments; the defaults are `atol = 1e-10`, `rtol = 0`. `evolve`, `evolve_bonddim` (including their methods taking `pp.PauliString` or `pp.PauliRotation` gates) and `evolve_deepcliffords` take the same keywords with the same defaults and pass them on to `apply!`.
  - Tests use exact Clifford angles, or angles at least 1e-3 away from every multiple of π/4. Angles with 1e-10 < |ϕ − mπ/4| ≲ 1e-8 are avoided in exact-regime tests: the new code takes them as generic, and the MPO apply then loses amplitudes of about 1e-9 (an inherited ITensorMPS precision floor, accepted by the user on 2026-11-12).
- **Circuit order.** `gates[1]` acts on the state first: |ψ_t⟩ = U_t ⋯ U_1 |ψ_0⟩. Angle vectors are in the same order as gate lists.
- **PauliPropagation.** `PauliRotation` with angle θ is exp(−iθP/2). `propagate` takes the gate list and the angle vector in the circuit order above and returns U†OU (checked numerically). A gate (P, ϕ) is passed to PP with θ = −2ϕ. A `PauliRotation` carries no sign, so (−P, ϕ) is passed as (P, θ = +2ϕ).
- **XXZ gate** on a bond (i, i+1): the three commuting gates (XX, ϕ_xy), (YY, ϕ_xy), (ZZ, ϕ_z), stored and applied in that order, all in the exp(iϕP) convention. `xxz_circuit(ϕ_xy, ϕ_z, t, N)` takes the angles in that order.
- **Bricklayer layer:** bonds (1,2), (3,4), … first, then (2,3), (4,5), …, each bond getting its full XX, YY, ZZ triple before the next bond. That is N−1 bonds and 3(N−1) gates per layer for any N (135 at N=46). Layer l ends at gate 3(N−1)·l. For other circuits the caller supplies the layer ends as gate indices.
- **Recording.** The `*_circuit_dynamics` functions return ⟨O⟩ at a list of record points.
  - Record points are the gate indices in `layer_ends`. With `layer_ends = nothing`, every gate is a record point.
  - Entry 1 of the result is the value at t = 0 (before any gate); entry l + 1 is the value after the gate at record point l.
  - PP evaluates the observable, and can stop on `Nmax`, only at record points.
  - `pauliprop_circuit_dynamics` continuing from a CAMPS state after `start_gate` gates returns only the record points after `start_gate`, with no t = 0 entry.
- **Layer-end disentangler** (`campssrc_circuit_dynamics`). At every record point, until the run is flagged as truncated, the CAMPS state is disentangled as in the old code: the `criterion` and `strategy` passed in, at most N sweeps (`max_iter = N`), stopping after the first sweep that gains ≤ `min_diff = 6`, with no SVD keywords. The truncation flag is set once norm² < (1 − thl) · norm²(t = 0); from then on the disentangler no longer runs.
  - **The flag is unreliable below thl = 1e-13** (user, 2026-11-16): rounding drift of the norm (about 1e-14) can set it with nothing truncated, so whether and when the disentangler runs then depends on floating-point details, and can differ between the old and new stacks. Behaviour at such thresholds is unspecified and not tested; tests that depend on the disentangler use thl ≥ 1e-11 and check that the norm drift stays far below thl.
- **Doping.** A mask is a Bool vector over the whole circuit, with (N−1)·n_layers entries; entry m refers to the m-th triple in gate order, i.e. gates 3m−2…3m. `xy_magic(phases, mask; magicphase)` replaces the XX and YY angles of every masked triple by `magicphase`; `z_magic` replaces the ZZ angle. `magicphase` is an angle ϕ in the exp(iϕP) convention; it has no fixed value (each script sets it). One call uses one mode for the whole mask.

## 3. Pauli strings and sums

- **PauliPropagation integer:** the term equals Σⱼ codeⱼ · 4^(j−1), with codes I=0, X=1, Y=2, Z=3; qubit j occupies bits 2(j−1) (low) and 2(j−1)+1 (high), for every key type.
  - **Key type** (`getinttype`): the smallest of `UInt8`, `UInt16`, `UInt32`, `UInt64`, `UInt128` that holds 2N bits (N ≤ 4, 8, 16, 32, 64); above 64 qubits `NTupleInteger{w}`, with w = ⌈2N/64⌉ 64-bit limbs, little-endian (limb 1 holds bits 0–63), so qubit j lies in limb ⌈j/32⌉. For example N = 10 → `UInt32`, N = 46 → `UInt128`, N = 70 → `NTupleInteger{3}`.
  - Every key type is an `Unsigned` integer with `<<`, `>>`, `&` and `Int(…)` of a value that fits. Tests take the type from the object, not from a literal, except in the pins of this rule. Coefficients are `Float64`.
  - Historical (PP 0.7.3): the key type was the unsigned integer with 8·⌈2N/8⌉ bits (N = 10 → `UInt24`, N = 46 → `UInt96`), the non-native widths created at runtime.
- **QuantumClifford:** (x, z) bits per qubit, X=(1,0), Y=(1,1), Z=(0,1). In terms of the PP bits, z = high and x = low ⊻ high.
- **Phases.** A QC `PauliOperator` with phase byte r stands for iʳ times the Hermitian string, so `P"Y"` is Y itself and a Hermitian operator has r ∈ {0, 2}. "The sign of a string" always means this sign.
- **Old converter** (`cmps.PauliSum(::pp.PauliSum)`; historical, deleted in Phase 6): one term per input term, coefficient unchanged, string with phase byte 0.
- **New converter** (`MagicTensors.QubitPauliSum(::pp.PauliSum)`): one term per input term, the same strings; coefficients converted to `ComplexF64` with zero imaginary part (input coefficients are real or complex). Each key is a `QubitPauli` wrapping a QC `PauliOperator` with phase byte 0, whose (x, z) bits follow the QC rule above, qubit j at position j, for any N (`UInt128` keys at N = 46 and `NTupleInteger` keys above 64 qubits included). No two keys share storage. Tests compare as dictionaries keyed by the length-N string of symbols (character j for qubit j), after dropping terms whose coefficient is exactly 0 on both sides; term order is unspecified.
- **Signed strings in a `QubitPauliSum`.** A term with key ±P and coefficient c stands for c·(±P): the sign of a key counts, wherever it is stored.

## 4. Clifford operators and the CAMPS frame

Throughout this section **C is the state-frame Clifford: |ψ⟩ = C |mps⟩.** All products are operator products.

- **QuantumClifford.** A `CliffordOperator` is defined by its tableau: the images C X_j C† and C Z_j C†. `C"…"` lists the X images first (X₁…X_N), then the Z images. With that definition:
  - `C * Stabilizer` and `apply!(stab, C)` give C P C†. QC has no `C * PauliOperator` (the vendored `CliffordMPS` added one, giving C P C†; it went with the old code in Phase 6). Tests write C P C† for a single string as `cliffconj(C, P)` = `(C * Stabilizer([P]))[1]` (test helper), and the code under test must not rely on `C * PauliOperator`;
  - `inv(C)` is the operator C† (= C⁻¹), so `cliffconj(inv(C), P)` and `apply!(stab, inv(C))` give C† P C;
  - `A * B` is the operator product A·B, so B acts first.

  Each of these is pinned by a direct hand-derived test (assumption pins). The pins show that the functions agree with the tableaus. That a tableau row means C X_j C† for the actual matrices of §1 (hence `tPhase` = diag(1, i), X → Y) is checked only by the dense reference, through physical expectation values.
- **`DisentangleCAMPS.apply(P, C)`** returns C†PC (the Heisenberg image).
- **Old frame** (historical). `CAMPS` stored `Cdag` = C⁻¹, so |ψ⟩ = Cdag⁻¹ |mps⟩ and ⟨ψ|P|ψ⟩ = ⟨mps| Cdag P Cdag⁻¹ |mps⟩. The old `paulinature` and `disentangler` took the state-frame C, i.e. `inv(ψ.Cdag)`, as their Clifford argument.
- **New frame** (checked against the pinned commit). `QubitCAMPS` holds |ψ⟩ = C |mps⟩. It stores either C (`daggered = false`) or C† (`daggered = true`, upstream's lazy flag); `get_clifford(ψ; daggered = d)` returns C for d = false and C† for d = true, whichever is stored. Tableau assertions are on the logical C (or C†), whatever the flag, and tests build non-trivial frames in both storages.
- **Frame update rules.** The state-frame column holds for `QubitCAMPS` (and held for the old stack); the old storage column (historical) is the same rule written for `Cdag`.

  | Operation on the state | State frame | Old storage |
  |---|---|---|
  | Clifford gate G | C ← G·C | Cdag ← Cdag·G⁻¹ |
  | Disentangled rotation, with (D, s) | C ← C·D, then exp(i·sϕ·X_{k+1}) on the MPS | Cdag ← D⁻¹·Cdag, then the same |
  | `:logical` rotation | C unchanged; exp(iϕ·C†PC) on the MPS | Cdag unchanged; the same |
  | `:trivial` rotation | nothing | nothing |
- **Frame-invariant tests** assert only on the physical state, through Pauli expectation values (and the stabilizer Rényi entropy, §5). **Frame-dependent tests** assert on `Cdag`/C, the MPS alone, k, or bond dimensions; they are tagged, and run on each stack with that stack's own rules (e.g. which angles are Clifford, §2). Since Phase 6 the only stack is `QubitCAMPS`.
- **Magic qubits.** With k magic qubits the MPS is |M⟩ = |m⟩₁…ₖ ⊗ |0⟩^(N−k), where |m⟩ is an arbitrary k-qubit state, in general entangled. k = 0 means |0…0⟩; k = N means no constraint.
  - k is tracked by the caller and returned by `apply!(ψ, k, P, ϕ)`. It increases by 1 only for a `:disentanglable` gate at an angle not recognised as Clifford; Clifford-angle gates go into the frame and leave k unchanged.
  - `addmagicstate!(ψ, k, α)`, with k the count before the call, multiplies qubit k+1 of the MPS by cos α · 1 + i sin α · X, i.e. turns its |0⟩ into exp(iαX)|0⟩. It does not change any counter. For a gate (P, ϕ) with disentangler sign s, α = s·ϕ.
- **`paulinature(k, C, P)`**, for Q = C†PC, ignoring the sign of Q: `:disentanglable` if Q has X or Y on a qubit > k; otherwise `:trivial` if Q is the identity on qubits 1…k; otherwise `:logical`. So at k = 0 every non-disentanglable Q is `:trivial`, and at k = N nothing is `:disentanglable`.
- **Disentangler.** `disentangler(k, C, P)` is tested at two levels. Let Q = C†PC have an X or Y on some qubit > k, and let s = ±1 be the sign of Q (§3). Its behaviour for any other Q is undefined; callers check `paulinature` first.
  - **Correctness properties** (frame-invariant; held on the old code and must hold on the new, on seeded random (k, C, P)). The returned (D, s) satisfy:
    1. D†QD = s · X_{k+1};
    2. D fixes every state |φ⟩₁…ₖ ⊗ |0…0⟩, up to one global phase. In tableau form, with G the group generated by Z_{k+1}, …, Z_N with sign +: D Z_j D† ∈ G for j > k, and D A D† ∈ A·G for A = X_j, Z_j with j ≤ k. (Requiring D A D† = A exactly is too strong: the pinned D fails it.)

    Then exp(iϕP)·C|M⟩ = (CD)·exp(i·sϕ·X_{k+1})|M⟩. Any D with these properties gives the right state.
  - **Pinned construction** (tagged implementation-specific). The properties leave D free, and the choice of D sets the frame, and with it later bond dimensions and the layer at which CAMPS-PP switches to PP. The construction is deterministic, so it is specified here and the tests pin the exact D:
    1. i = the smallest index > k where Q has X or Y.
    2. If Q_i = Y: D₁ = S_i, the phase gate diag(1, i) on qubit i (S†YS = X). Otherwise D₁ = 1.
    3. If i ≠ k+1: D₂ = SWAP(i, k+1). Otherwise D₂ = 1.
    4. Let Q′ = (D₁D₂)† Q (D₁D₂), which has X on qubit k+1 and the same sign as Q. D₃ = the product, over all j ≠ k+1 (both j ≤ k and j > k+1) where Q′ is not the identity, of controlled-σ_j = |0⟩⟨0|_{k+1} ⊗ 1 + |1⟩⟨1|_{k+1} ⊗ σ_j, with σ_j ∈ {X, Y, Z} the unsigned Pauli of Q′ on qubit j (CNOT, CY, CZ; control k+1, target j). These commute.
    5. D = D₁ D₂ D₃, so D₃ acts first on a state.

    Expected tableaus are worked out by hand from steps 1–5 and must also pass properties 1–2. `build_D`, `swapqubits` and `reducetoX` are pinned only through the full D; the intermediate sign returned by `reducetoX` is bookkeeping and not part of the spec.
  - **Port.** If the reimplementation keeps this construction, the pins carry over unchanged. If it uses the upstream choice of D, the pins are rewritten deliberately and the resource changes are reviewed in the old-vs-new comparison.

## 5. Observables and entropies

- **`expectation(ψ, O)`** returns ⟨ψ|O|ψ⟩ for the physical state |ψ⟩ = C|mps⟩, i.e. ⟨mps|C†OC|mps⟩, not divided by ⟨mps|mps⟩. O is a signed Pauli string or a Pauli sum. Tests compare the real part and require |imaginary part| ≤ 1e-10.
- **The pp seam** `expectation(ψ::QubitCAMPS, O::pp.PauliSum; method = :perterm, normalized = false, verbose = false)` (post-port item 4; decided by the user 2026-11-19).
  - It returns Σ c_P ⟨mps|C†PC|mps⟩ over the terms (P, c_P) of O, i.e. ⟨ψ|O|ψ⟩ as above, as a `ComplexF64`. With `normalized = true` it returns that value divided by ⟨mps|mps⟩ (= ⟨ψ|ψ⟩). The default is `false`.
  - Coefficients may be real or complex. With complex coefficients the value is complex, and tests compare the whole complex number. The identity term contributes c·⟨mps|mps⟩. An empty sum gives 0.
  - `method = :perterm` (the default) evaluates the sum term by term. `method = :mpo` is the earlier route (conjugate → OpSum → MPO → `inner`), kept for checks. The two agree to 1e-12 on the tested sums (the exact-regime cases and the T8 hand-off sums). Any other `method` raises `ArgumentError`. How `:perterm` orders the terms or reuses environments is not part of the spec.
  - It accepts an MPS in any gauge (any orthogonality centre, or none) and either frame storage (§4).
  - ψ is unchanged: after the call, the MPS tensors (with their site and link indices) and the stored frame, `daggered` flag included, are bitwise equal to what they were before.
  - It runs serially, with no threads inside one call, so repeated calls on the same input return bitwise-equal values.
  - `verbose = true` prints progress and returns the same value.
- **Transferred magnetization.** `transferredmagnetization(N, init_ones)`, with L = {1, …, ⌊N/2⌋}, R = the other qubits, n₁ᴸ(0) = |init_ones ∩ L| and n₀ᴿ(0) = |R ∖ init_ones|, returns the `pp.PauliSum`

  O = −½ Σ_{j∈L} Z_j + ½ Σ_{j∈R} Z_j + (N/2 − n₁ᴸ(0) − n₀ᴿ(0)) · 1.

  - This is (ones gained by L) + (zeros gained by R), and ⟨O⟩ = 0 on the initial state. The identity term is always present, also when its coefficient is 0.
  - For a circuit that conserves ΣZ, ⟨O⟩ = 2 × (net number of 1s moved from R to L): one excitation crossing right → left gives +2, left → right gives −2. A single flip that breaks the conservation gives ±1.
  - The factor 2 is intended (user, 2026-11-11): it is the mathematical definition in Rosenberg et al.
- **Entanglement entropy:** von Neumann, log base 2: S = −Σ p log₂ p over the squared Schmidt values p of a bond, with 0·log 0 = 0. It is defined for a normalised state (Σp = 1); the functions do not normalise, and tests use normalised states. `entanglement_entropy` returns N−1 values, element b for the bond between sites b and b+1. On a CAMPS it is the entropy of the MPS part, which depends on the frame, so it is a frame-dependent quantity. The code has no function for the entropy of the physical state.
- **Stabilizer Rényi entropy:** for a normalised state, with Ξ_P = ⟨P⟩²/2ᴺ over all 4ᴺ unsigned strings (ΣΞ = 1): M_α = (1−α)⁻¹ log₂ Σ_P Ξ_Pᵅ − N; α = 1 is −Σ Ξ log₂ Ξ − N; the linear entropy is 1 − 2ᴺ Σ Ξ². Not divided by N (scripts do that themselves). M_α = 0 for stabilizer states; the single-qubit state exp(iπ/8·X)|0⟩ has M₂ = log₂(4/3). It is invariant under Cliffords, so the value for the MPS part equals that of the physical state, and it is frame-invariant.
- **Entropy functions:** MagicTensors' `stabilizer_entropy` (exact, and sampled with a sample count) and `entanglement_entropy`, tested against these definitions; both take a plain MPS. The old `sEntropy`/`eEntropys!` (CliffordMPS) followed the same definitions and were tested the same way until Phase 6.

## 6. Tolerances

All tolerances are absolute.

- **Exact regime** (no truncation): 1e-10 on each expectation value, whether of a single string or of a sum; 1e-10 on entropies.
- **Old vs new,** on replayed circuits. A *sample* is one circuit realisation (initial bitstring plus doping mask); the quantity is ⟨O⟩ at the end of each layer.
  - CAMPS-PP layers that come before the switch to PP in both the old and the new run: |Δ⟨O⟩| ≤ 1e-10.
  - Every other layer and method: |Δ⟨O⟩| ≤ δ = 1e-6, per sample and per layer. A sample above δ is investigated individually; δ is not loosened without the user's written approval.
- **Norm drift in the truncating CAMPS method.** `apply!(ψ, P, ϕ, thl)` sends every non-Clifford gate through an SVD and never renormalises, and `expectation` is not normalised (§5), so ⟨O⟩ carries the rounding drift of the norm (measured: norm² off by 4e-12 after 310 gates at N = 46, i.e. 1.3e-10 on ⟨O⟩ ≈ 30). For this method, over long circuits, tests require |norm² − 1| ≤ 1e-10 and compare ⟨O⟩ up to that drift, or after dividing by norm².
- **Renormalisation in `mps_circuit_dynamics`.** Unlike the truncating CAMPS method, it truncates after every gate at `cutoff = thl` and then renormalises, as the old code did, so the state it returns has norm 1 and the ⟨O⟩ it records are those of a normalised state. Tests require |norm − 1| ≤ 1e-12 on the returned state.
- **Truncated regime.** No absolute error bound is specified. Tests require only that the error vanishes without truncation and, on the configuration tested, that it does not grow as the threshold is tightened. The latter is not a general property of the methods (for the MPS method, one small configuration gives an error of 0.016 at thl = 0.1 and 0.028 at thl = 0.03), so a threshold-ladder test uses a configuration where it has been checked to hold.
- Per-method parameters stay separate (χ_camps, thl_pp, Nmax_pp).

## Decisions (all approved by the user, 2026-11-11)

1. **Rotation convention** (§2): exp(iϕP) for all signed ϕ. Catching Clifford angles in dedicated branches is fine and stays; only the two bug classes listed in §2 are wrong. The branches for 0, π, +π/4 and −3π/4 agree with the convention.
2. **Transferred magnetization** (§5): the operator with the factor 2 is the intended one; only the docstring is wrong.
3. **Entropy normalisation** (§5): log base 2 for both entropies, M_α as written, no division by N inside the functions.
4. **Disentangler** (§4): properties 1–2 plus the pinned construction, with the swap in step 3.
5. **New frame** (§4): provisional until the upstream commit is pinned. **Resolved 2026-11-13:** checked against the pinned commit and no longer provisional.

## Change log

- 2026-11-11: drafted; decisions 1–5 approved.
- 2026-11-11: older tests that expected CPC† from `apply(P, C)`, and one that expected a disentangler without the swap, were stale and have been replaced.
- 2026-11-11: revised after two independent reviews (one as a test author, one as the port implementer; each read only this page). Added: the terms list and one symbol per meaning; the state-frame definition of C and the frame update table; the exact scope of the two bug classes; mask indexing and `magicphase`; the PP key-type rule and the handling of signed gates; the tableau form of property 2; the definition of controlled-σ_j and of the sign s; magic-qubit bookkeeping; `expectation` on the physical state; the transferred-magnetization formula with ⌊N/2⌋ and the identity term; entropy normalisation and bond order; absolute tolerances and the definition of a sample. Removed "full-state entropies" from the frame-invariant list (the code has no such function).
- 2026-11-11: after the adversarial test review, bug class 2 in §2 was narrowed from "the frame is not a Pauli operator" to "C†PC ≠ ±P" (the reviewer showed a non-Pauli frame where ±π/2 is correct).
- 2026-11-11: added, with the user's approval: the recording conventions (§2).
- 2026-11-11: added, with the user's approval: the domain-wall ensemble and its orientation (§1); the norm-drift allowance for the truncating CAMPS method and the absence of a bound for the truncated regime (§6).
- 2026-11-13, with the user's approval: updated for the port. §2: the old-type-only scope of the known bugs; the new code's Clifford-angle rule and its defaults (atol 1e-10, rtol 0, every multiple of π/4; decided 2026-11-12); angles to avoid because of the ~1e-9 MPO precision floor. §3: the new converter's key layout, no longer provisional; signed keys. §4: the new frame, checked against the pinned commit, with both storages; frame-dependent tests run per stack. §5: the two entropy function families. Terms: the pinned MagicTensors commit; QuantumClifford is 0.11.4 (the Manifest version), not 0.11.3. No convention changed.
- 2026-11-16, with the user's approval, after the T9 review: §2: invalid input (±i strings, non-finite angles) rejected with `ArgumentError`; the `atol`/`rtol` keywords also cover the truncating `apply!` and are passed on by `evolve`, `evolve_bonddim` and `evolve_deepcliffords`; the layer-end disentangler of `campssrc_circuit_dynamics` (old-code parameters) and its truncation flag, unreliable below thl = 1e-13. §6: `mps_circuit_dynamics` renormalises; the truncated-regime monotonicity holds per tested configuration, not in general. No convention changed.
- 2026-11-17, with the user's approval, after Phase 6 deleted the old code: old-code items marked historical (terms, §2 known bugs, §3 old converter, §4 old frame and old storage column, §5 old entropy functions); §4: QC has no `C * PauliOperator` (it was the vendored CliffordMPS's), tests use the helper `cliffconj`; §5: `entanglement_entropy` named for the bond order. No convention changed.
- 2026-11-19, with the user's approval: PP ≥ 0.9 key types (§3, Terms). The bit layout is unchanged; no convention changed.
- 2026-11-19, with the user's approval: §5, the pp seam for post-port item 4 (per-term evaluation; the MPO route kept as `method = :mpo`; the `normalized` keyword; ψ unchanged; serial). Values are unchanged up to rounding; no convention changed.
- Open follow-ups outside this page: none. (Done: the `transferredmagnetization` docstring, Phase 2; the key layout, the new-frame check and the Clifford-angle defaults, 2026-11-13.)
