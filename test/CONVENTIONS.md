# Convention spec (T2)

**Status: approved by the user on 2026-11-11** (five decisions, listed at the end). The new-frame items in §4 are approved as provisional. Revised the same day after two independent ambiguity reviews; the revision changed wording and added definitions, not conventions (see the change log).

Every expected value in the tests comes from this page, the dense reference or the hand-derived golden set, never from the implementation.

**Terms used on this page**
- *Old code*: the current packages (`DisentangleCAMPS`, `CampsPP`, on the vendored `CliffordMPS`). *New code*: the port to MagicTensors.jl.
- *Dense reference*: `DisentangleCAMPS/test/DenseReference.jl`, an exact state-vector simulator written without reading the implementation.
- *Golden set / assumption pins*: hand-derived checks in `DisentangleCAMPS/test/TestGoldenSet.jl` and `TestEvolve.jl`.
- *PP*: PauliPropagation.jl (0.7.3). *QC*: QuantumClifford.jl (0.11.3).
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
- **Known bugs in the old `apply!`** (both methods). Only these cases are `@test_broken` on the old code, and they must pass after the port:
  1. ϕ = −π/4 and ϕ = +3π/4 are wrong in every frame, including the identity frame;
  2. ϕ = ±π/2 is wrong when the frame is not a Pauli operator (the identity counts as Pauli), i.e. after a disentangler or a π/4 gate. With a Pauli frame it is correct.
- **Clifford angles and tolerance.** The physical state must equal exp(iϕP)|ψ⟩ for every real ϕ. Which angles an implementation recognises as Clifford affects only frame-dependent quantities (k, the frame, bond dimensions).
  - The old code recognises |ϕ| ∈ {0, π/4, π/2, 3π/4, π}; other multiples of π/4 (5π/4, 2π, …) take the generic branch. It uses `atol = 1e-6` at 0 and the default `≈` relative tolerance elsewhere.
  - The port makes this uniform: optional `atol` and `rtol` parameters applying to every Clifford angle (user, 2026-11-11); defaults are chosen then.
  - Until then tests use either exact Clifford angles or angles at least 1e-3 away from every multiple of π/4.
- **Circuit order.** `gates[1]` acts on the state first: |ψ_t⟩ = U_t ⋯ U_1 |ψ_0⟩. Angle vectors are in the same order as gate lists.
- **PauliPropagation.** `PauliRotation` with angle θ is exp(−iθP/2). `propagate` takes the gate list and the angle vector in the circuit order above and returns U†OU (checked numerically). A gate (P, ϕ) is passed to PP with θ = −2ϕ. A `PauliRotation` carries no sign, so (−P, ϕ) is passed as (P, θ = +2ϕ).
- **XXZ gate** on a bond (i, i+1): the three commuting gates (XX, ϕ_xy), (YY, ϕ_xy), (ZZ, ϕ_z), stored and applied in that order, all in the exp(iϕP) convention. `xxz_circuit(ϕ_xy, ϕ_z, t, N)` takes the angles in that order.
- **Bricklayer layer:** bonds (1,2), (3,4), … first, then (2,3), (4,5), …, each bond getting its full XX, YY, ZZ triple before the next bond. That is N−1 bonds and 3(N−1) gates per layer for any N (135 at N=46). Layer l ends at gate 3(N−1)·l. For other circuits the caller supplies the layer ends as gate indices.
- **Doping.** A mask is a Bool vector over the whole circuit, with (N−1)·n_layers entries; entry m refers to the m-th triple in gate order, i.e. gates 3m−2…3m. `xy_magic(phases, mask; magicphase)` replaces the XX and YY angles of every masked triple by `magicphase`; `z_magic` replaces the ZZ angle. `magicphase` is an angle ϕ in the exp(iϕP) convention; it has no fixed value (each script sets it). One call uses one mode for the whole mask.

## 3. Pauli strings and sums

- **PauliPropagation integer:** the term equals Σⱼ codeⱼ · 4^(j−1), with codes I=0, X=1, Y=2, Z=3; qubit j occupies bits 2(j−1) (low) and 2(j−1)+1 (high). The key type is the unsigned integer with 8·⌈2N/8⌉ bits (N ≤ 4 → `UInt8`, N = 10 → `UInt24`, N = 46 → `UInt96`); tests take it from the object, not from a literal. Coefficients are `Float64`.
- **QuantumClifford:** (x, z) bits per qubit, X=(1,0), Y=(1,1), Z=(0,1). In terms of the PP bits, z = high and x = low ⊻ high.
- **Phases.** A QC `PauliOperator` with phase byte r stands for iʳ times the Hermitian string, so `P"Y"` is Y itself and a Hermitian operator has r ∈ {0, 2}. "The sign of a string" always means this sign.
- **Old converter** (`cmps.PauliSum(::pp.PauliSum)`): one term per input term, coefficient unchanged, string with phase byte 0.
- **New converter (provisional until the upstream commit is pinned):** `pp.PauliSum` → upstream `QubitPauliSum`. Same strings; coefficients converted to `ComplexF64` with zero imaginary part; input coefficients are real. Tests compare as dictionaries keyed by the length-N tuple of symbols (entry j for qubit j), after dropping terms whose coefficient is exactly 0 on both sides; term order is unspecified. The upstream key layout is fixed in this spec once the commit is pinned.

## 4. Clifford operators and the CAMPS frame

Throughout this section **C is the state-frame Clifford: |ψ⟩ = C |mps⟩.** All products are operator products.

- **QuantumClifford.** A `CliffordOperator` is defined by its tableau: the images C X_j C† and C Z_j C†. `C"…"` lists the X images first (X₁…X_N), then the Z images. With that definition:
  - `C * Stabilizer`, `C * PauliOperator` (a method added by the vendored `CliffordMPS`) and `apply!(stab, C)` all give C P C†;
  - `inv(C)` is the operator C† (= C⁻¹), so `inv(C) * P` and `apply!(stab, inv(C))` give C† P C;
  - `A * B` is the operator product A·B, so B acts first.

  Each of these is pinned by a direct hand-derived test (assumption pins). The pins show that the functions agree with the tableaus. That a tableau row means C X_j C† for the actual matrices of §1 (hence `tPhase` = diag(1, i), X → Y) is checked only by the dense reference, through physical expectation values.
- **`DisentangleCAMPS.apply(P, C)`** returns C†PC (the Heisenberg image).
- **Old frame.** `CAMPS` stores `Cdag` = C⁻¹, so |ψ⟩ = Cdag⁻¹ |mps⟩ and ⟨ψ|P|ψ⟩ = ⟨mps| Cdag P Cdag⁻¹ |mps⟩. The old `paulinature` and `disentangler` take the state-frame C, i.e. `inv(ψ.Cdag)`, as their Clifford argument.
- **New frame (provisional; approved as such).** Upstream `QubitCAMPS` stores C with |ψ⟩ = C |mps⟩, plus a lazy `daggered` flag. This is to be re-checked against the pinned upstream commit before any test relies on it. Tableau assertions are on the logical C, whatever the flag.
- **Frame update rules.**

  | Operation on the state | State frame | Old storage |
  |---|---|---|
  | Clifford gate G | C ← G·C | Cdag ← Cdag·G⁻¹ |
  | Disentangled rotation, with (D, s) | C ← C·D, then exp(i·sϕ·X_{k+1}) on the MPS | Cdag ← D⁻¹·Cdag, then the same |
  | `:logical` rotation | C unchanged; exp(iϕ·C†PC) on the MPS | Cdag unchanged; the same |
  | `:trivial` rotation | nothing | nothing |
- **Frame-invariant tests** assert only on the physical state, through Pauli expectation values (and the stabilizer Rényi entropy, §5). **Frame-dependent tests** assert on `Cdag`/C, the MPS alone, k, or bond dimensions; they are tagged and rewritten for the new code.
- **Magic qubits.** With k magic qubits the MPS is |M⟩ = |m⟩₁…ₖ ⊗ |0⟩^(N−k), where |m⟩ is an arbitrary k-qubit state, in general entangled. k = 0 means |0…0⟩; k = N means no constraint.
  - k is tracked by the caller and returned by `apply!(ψ, k, P, ϕ)`. It increases by 1 only for a `:disentanglable` gate at an angle not recognised as Clifford; Clifford-angle gates go into the frame and leave k unchanged.
  - `addmagicstate!(ψ, k, α)`, with k the count before the call, multiplies qubit k+1 of the MPS by cos α · 1 + i sin α · X, i.e. turns its |0⟩ into exp(iαX)|0⟩. It does not change any counter. For a gate (P, ϕ) with disentangler sign s, α = s·ϕ.
- **`paulinature(k, C, P)`**, for Q = C†PC, ignoring the sign of Q: `:disentanglable` if Q has X or Y on a qubit > k; otherwise `:trivial` if Q is the identity on qubits 1…k; otherwise `:logical`. So at k = 0 every non-disentanglable Q is `:trivial`, and at k = N nothing is `:disentanglable`.
- **Disentangler.** `disentangler(k, C, P)` is tested at two levels. Let Q = C†PC have an X or Y on some qubit > k, and let s = ±1 be the sign of Q (§3). Its behaviour for any other Q is undefined; callers check `paulinature` first.
  - **Correctness properties** (frame-invariant; must hold on the old and the new code, on seeded random (k, C, P)). The returned (D, s) satisfy:
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
- **Transferred magnetization.** `transferredmagnetization(N, init_ones)`, with L = {1, …, ⌊N/2⌋}, R = the other qubits, n₁ᴸ(0) = |init_ones ∩ L| and n₀ᴿ(0) = |R ∖ init_ones|, returns the `pp.PauliSum`

  O = −½ Σ_{j∈L} Z_j + ½ Σ_{j∈R} Z_j + (N/2 − n₁ᴸ(0) − n₀ᴿ(0)) · 1.

  - This is (ones gained by L) + (zeros gained by R), and ⟨O⟩ = 0 on the initial state. The identity term is always present, also when its coefficient is 0.
  - For a circuit that conserves ΣZ, ⟨O⟩ = 2 × (net number of 1s moved from R to L): one excitation crossing right → left gives +2, left → right gives −2. A single flip that breaks the conservation gives ±1.
  - The factor 2 is intended (user, 2026-11-11): it is the mathematical definition in Rosenberg et al.
- **Entanglement entropy:** von Neumann, log base 2: S = −Σ p log₂ p over the squared Schmidt values p of a bond, with 0·log 0 = 0. It is defined for a normalised state (Σp = 1); the functions do not normalise, and tests use normalised states. `eEntropys!` returns N−1 values, element b for the bond between sites b and b+1. On a CAMPS it is the entropy of the MPS part, which depends on the frame, so it is a frame-dependent quantity. The code has no function for the entropy of the physical state.
- **Stabilizer Rényi entropy:** for a normalised state, with Ξ_P = ⟨P⟩²/2ᴺ over all 4ᴺ unsigned strings (ΣΞ = 1): M_α = (1−α)⁻¹ log₂ Σ_P Ξ_Pᵅ − N; α = 1 is −Σ Ξ log₂ Ξ − N; the linear entropy is 1 − 2ᴺ Σ Ξ². Not divided by N (scripts do that themselves). M_α = 0 for stabilizer states; the single-qubit state exp(iπ/8·X)|0⟩ has M₂ = log₂(4/3). It is invariant under Cliffords, so the value for the MPS part equals that of the physical state, and it is frame-invariant.

## 6. Tolerances

All tolerances are absolute.

- **Exact regime** (no truncation): 1e-10 on each expectation value, whether of a single string or of a sum; 1e-10 on entropies.
- **Old vs new,** on replayed circuits. A *sample* is one circuit realisation (initial bitstring plus doping mask); the quantity is ⟨O⟩ at the end of each layer.
  - CAMPS-PP layers that come before the switch to PP in both the old and the new run: |Δ⟨O⟩| ≤ 1e-10.
  - Every other layer and method: |Δ⟨O⟩| ≤ δ = 1e-6, per sample and per layer. A sample above δ is investigated individually; δ is not loosened without the user's written approval.
- **Norm drift in the truncating CAMPS method.** `apply!(ψ, P, ϕ, thl)` sends every non-Clifford gate through an SVD and never renormalises, and `expectation` is not normalised (§5), so ⟨O⟩ carries the rounding drift of the norm (measured: norm² off by 4e-12 after 310 gates at N = 46, i.e. 1.3e-10 on ⟨O⟩ ≈ 30). For this method, over long circuits, tests require |norm² − 1| ≤ 1e-10 and compare ⟨O⟩ up to that drift, or after dividing by norm².
- **Truncated regime.** No absolute error bound is specified. Tests require only that the error does not grow as the threshold is tightened and that it vanishes without truncation.
- Per-method parameters stay separate (χ_camps, thl_pp, Nmax_pp).

## Decisions (all approved by the user, 2026-11-11)

1. **Rotation convention** (§2): exp(iϕP) for all signed ϕ. Catching Clifford angles in dedicated branches is fine and stays; only the two bug classes listed in §2 are wrong. The branches for 0, π, +π/4 and −3π/4 agree with the convention.
2. **Transferred magnetization** (§5): the operator with the factor 2 is the intended one; only the docstring is wrong.
3. **Entropy normalisation** (§5): log base 2 for both entropies, M_α as written, no division by N inside the functions.
4. **Disentangler** (§4): properties 1–2 plus the pinned construction, with the swap in step 3.
5. **New frame** (§4): provisional until the upstream commit is pinned.

## Change log

- 2026-11-11: drafted; decisions 1–5 approved.
- 2026-11-11: older tests that expected CPC† from `apply(P, C)`, and one that expected a disentangler without the swap, were stale and have been replaced.
- 2026-11-11: revised after two independent reviews (one as a test author, one as the port implementer; each read only this page). Added: the terms list and one symbol per meaning; the state-frame definition of C and the frame update table; the exact scope of the two bug classes; mask indexing and `magicphase`; the PP key-type rule and the handling of signed gates; the tableau form of property 2; the definition of controlled-σ_j and of the sign s; magic-qubit bookkeeping; `expectation` on the physical state; the transferred-magnetization formula with ⌊N/2⌋ and the identity term; entropy normalisation and bond order; absolute tolerances and the definition of a sample. Removed "full-state entropies" from the frame-invariant list (the code has no such function).
- 2026-11-11: added, with the user's approval: the domain-wall ensemble and its orientation (§1); the norm-drift allowance for the truncating CAMPS method and the absence of a bound for the truncated regime (§6).
- Open follow-ups outside this page: the `transferredmagnetization` docstring (porting session); the upstream `QubitPauliSum` key layout and the new-frame check (when the commit is pinned); default `atol`/`rtol` for Clifford angles (port).
