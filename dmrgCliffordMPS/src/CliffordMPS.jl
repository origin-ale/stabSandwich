"""
CliffordMPS provides Clifford-augmented tensor-network tools for qubit systems.

The package combines `ITensors`, `ITensorMPS`, and `QuantumClifford` to work with:

- dressed MPS states stored as `CAMPS`,
- sparse Pauli operators and Pauli-sum Hamiltonians,
- Clifford gate indexing, circuit synthesis, and conjugation,
- entanglement and stabilizer-magic diagnostics,
- disentangling heuristics and local Clifford rotations,
- DMRG and TDVP on the non-Clifford tensor-network part.

The exported API is organized by feature rather than by implementation file. See the package README for a full user-level overview.
"""
module CliffordMPS

export density_matrix
export sigmas
include("ManipulateStab.jl")

export applyGate!
export applySwaps!
export spectrum!
export spectra!
export eEntropy
export eEntropy!
export eEntropys!
export random_qubit_mps
include("ManipulateMPS.jl")

export sEntropy
include("StabilizerEntropy.jl")

export clifford_gate
export clifford_number
export Clifford1Index
export Clifford1IndexSet
export Clifford2Index
export Clifford2IndexSet
export random_clifford_circuit
include("CliffordGate.jl")

export PauliSum
export expectation
export updateMPO!
include("PauliSum.jl")

export CAMPS
export transform!
export random_camps
export dmrg!
export tdvp!
include("ClifAugMPS.jl")

export DisentangleCriterion
export DisentangleStrategy
export disentangler
export disentangle!
export local_rotation!
include("Disentangler.jl")

export QOp
export applyQOp!
include("QuantumOperations.jl")

end

