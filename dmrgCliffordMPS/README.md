>[!NOTE]
> Copilot-generated documentation, take with a grain of salt.

# CliffordMPS

CliffordMPS is a Julia package for working with quantum states that are represented as an MPS dressed by a Clifford frame. It combines `ITensors`, `ITensorMPS`, and `QuantumClifford` so that Clifford circuits, Pauli Hamiltonians, stabilizer observables, tensor-network evolution, and disentangling heuristics can be used in one workflow.

The central representation is `CAMPS`, short for Clifford-Augmented MPS. A `CAMPS` value stores an MPS together with a Clifford operator `Cdag`, and the physical state is interpreted as the dressed state obtained from that pair. In practice, the package lets you move Clifford gates between the frame and the tensor network, measure local observables, estimate stabilizer magic, and run DMRG/TDVP on the non-Clifford remainder.

## Installation

Activate the package project and instantiate it in Julia:

```julia
using Pkg
Pkg.activate("/Users/alessandro/Programmi/magiCAMPS/dmrgCliffordMPS")
Pkg.instantiate()
```

The project was developed against Julia 1.10.x. Newer compatible Julia releases should also work, subject to the compatibility of the package dependencies.

## Quick Start

```julia
using CliffordMPS
using QuantumClifford

N = 12
camps = random_camps(N)

# Expectation values on the dressed state.
ex = expectation(camps, P"Z", [1])

# Apply a local Clifford gate.
applyQOp!(camps, QOp(:H), 1)

# Apply a random Clifford circuit and then search for disentangling moves.
applyQOp!(camps, QOp(:randCliffCircuit))
diff, schedule, cost = disentangle!(camps)

# One-qubit rotation that aligns each local Bloch vector with the computational basis.
distances = local_rotation!(camps)

# Stabilizer magic entropy of the MPS part.
magic = sEntropy(camps)
```

## Package Model

### `CAMPS`

`CAMPS` couples two objects:

- an `ITensorMPS.MPS` storing the non-Clifford tensor network part,
- a `QuantumClifford.CliffordOperator` storing the Clifford dressing.

The package provides methods that apply gates either to the full dressed state or by pushing Clifford operations into the frame when that is cheaper.

### `PauliSum`

`PauliSum` is a sparse Hamiltonian-like container built from a `Stabilizer` basis and a coefficient vector. It can be used to:

- compute expectation values on `Stabilizer` states and on `CAMPS`,
- convert to an `ITensorMPS.OpSum` or `MPO`,
- apply Clifford conjugation to the operator itself.

### Entropy and magic

The package includes two families of diagnostics:

- ordinary bipartite entanglement entropy from MPS singular values,
- stabilizer Renyi entropies and linear entropy based on Pauli expectations or sampling.

The stabilizer-entropy routines follow the constructions of Leone, Oliviero, and Hamma for exact evaluation and Lami and Collura for sampling-based estimation.

## Main API

### State construction and manipulation

- `CAMPS(mps)`, `CAMPS(Cdag)`, `CAMPS(N)` construct dressed states.
- `random_camps(N; debth, linkdims, type)` creates a random dressed state.
- `transform!(camps, ind, sites...)` applies a Clifford gate to both components while keeping the representation consistent.
- `applyGate!(camps, h::PauliSum; svd_kwargs...)` applies a Pauli-sum operator to the state.
- `applyGate!(camps, ind::Clifford1Index, n)` and `applyGate!(camps, ind::Clifford2Index, nL, nR)` apply indexed Clifford gates.

### Expectation values and density matrices

- `expectation(camps, h::PauliSum)` evaluates a Pauli Hamiltonian on a dressed state.
- `expectation(camps, p::PauliOperator[, sites])` evaluates a Pauli string.
- `density_matrix(camps, sites)` and `density_matrix(psi::Stabilizer, sites)` return reduced density matrices.

### Tensor-network algorithms

- `dmrg!(camps, h::PauliSum; dmrg_kwargs...)` runs DMRG on the dressed Hamiltonian.
- `tdvp!(camps, h::PauliSum, t; tdvp_kwargs...)` runs TDVP time evolution.
- `random_qubit_mps(N; linkdims, type)` builds a random qubit MPS for experiments.
- `applyGate!`, `applySwaps!`, `spectrum!`, `spectra!`, `eEntropy!`, and `eEntropys!` are the lower-level MPS utilities used by the algorithms above.

### Clifford circuit utilities

- `Clifford1Index`, `Clifford2Index` and their `...Set` variants provide a structured numbering scheme for Clifford gates.
- `clifford_gate(ind)` returns the corresponding `QuantumClifford.CliffordOperator`.
- `random_clifford_circuit(N, depth; clifford_set)` builds random Clifford circuits on `N` sites.

### Disentangling

- `DisentangleCriterion(:entangle)` and `DisentangleCriterion(:chi3)` select the cost function.
- `DisentangleStrategy(:full)`, `:radius`, `:brickwork`, `:snake`, and `:none` choose the scan pattern.
- `disentangler(mps, nL, nR; ...)` searches for the best two-site Clifford move.
- `disentangle!(camps, strategy; ...)` applies the selected schedule to a `CAMPS` state.
- `local_rotation!(camps)` applies one-qubit Cliffords that align local Bloch vectors with the computational axis.

### Quantum operations

- `QOp(x)` constructs a typed quantum operation token.
- `applyQOp!(camps, qop, ...)` is the high-level interface for gates, projective measurements, normalization, scaling, and random circuit layers.

## Examples of Common Workflows

### 1. Build a Hamiltonian from Pauli terms

```julia
using CliffordMPS
using QuantumClifford

N = 8
h = PauliSum(N)
for n in 1:N-1
	push!(h, (1.0, PauliOperator(P"XX", N, [n, n+1])))
end
```

### 2. Evolve and disentangle a dressed state

```julia
camps = random_camps(10)
energy = dmrg!(camps, h; nsweeps=2, cutoff=1e-8)
diff, schedule, cost = disentangle!(camps, DisentangleStrategy(:snake), 5)
```

### 3. Estimate stabilizer magic

```julia
magic_exact = sEntropy(camps, α=2.0)
magic_sampled = sEntropy(camps.mps, 128; α=2.0)
```

## Notes

- `random_qubit_mps(N; type="ITensors")` returns a normalized MPS compatible with ITensors' site conventions.
- Several routines expect the MPS to be orthogonally centered before applying local measurements or Pauli-sampling estimates. The helper methods do this internally when needed.
- `VectorizedMPDO` is included for an auxiliary test path and requires an even-length underlying MPS.

## References

- Leone, Oliviero, Hamma, Phys. Rev. Lett. 128, 050402 (2022). Stabilizer Renyi entropy from Pauli expectation values.
- Lami, Collura, arXiv:2303.05536. Sampling-based estimators for stabilizer entropies.

## Package Contents

The main implementation lives under `src/`:

- `CliffordMPS.jl` package entry point and exports.
- `ClifAugMPS.jl` dressed states and DMRG/TDVP helpers.
- `CliffordGate.jl` Clifford indexing and gate construction.
- `Disentangler.jl` disentangling heuristics.
- `ManipulateMPS.jl` MPS gate application and entropy utilities.
- `ManipulateStab.jl` stabilizer utilities and density matrices.
- `PauliSum.jl` Pauli Hamiltonian representation.
- `QuantumOperations.jl` high-level operation tokens.
- `StabilizerEntropy.jl` stabilizer entropy estimators.
- `VectorizedMPO.jl` auxiliary vectorized MPDO container.

