# Shared helpers for the DisentangleCAMPS tests (included once from runtests.jl).

using QuantumClifford
using Test
import ITensorMPS

"""
C P C† for a CliffordOperator C and a Pauli string P, sign included. QuantumClifford has no
`C * P` (the old vendored CliffordMPS defined it; pinned absent in TestGoldenSet.jl), so this
goes through QC's `C * Stabilizer`, pinned as C P C† in TestGoldenSet.jl.
"""
cliffconj(C::CliffordOperator, P::PauliOperator) = (C * Stabilizer([P]))[1]

"""
The gate `gate` acting on qubits `inds` of an N-qubit register, as a CliffordOperator: local
qubit a of the gate is qubit inds[a]. Built from the gate's tableau (images of X_a, Z_a, via
`cliffconj`) placed with QuantumClifford's `embed(N, inds, P)`; generators outside `inds` map
to themselves. Contract pinned in TestDisentangler.jl ("assumption pins: controlled gates and
embedding").
"""
function embed(gate::CliffordOperator, N, inds)
  n = nqubits(gate)
  length(inds) == n && allunique(inds) && all(1 .<= inds .<= N) ||
    throw(ArgumentError("need $n distinct qubits in 1:$N, got $inds"))
  image(local_gen, j) = (a = findfirst(==(j), inds);
    isnothing(a) ? local_gen(N, j) : QuantumClifford.embed(N, inds, cliffconj(gate, local_gen(n, a))))
  rows = vcat([image(single_x, j) for j in 1:N], [image(single_z, j) for j in 1:N])
  return CliffordOperator(QuantumClifford.Tableau(rows))
end

"""
Nature of a rotation about Q = C†PC with k magic qubits, straight from CONVENTIONS.md §4
(the sign of Q is ignored). Independent of `DisentangleCAMPS.paulinature`.
"""
function spec_nature(k, Q)
  any(xbit(Q)[k+1:end]) && return :disentanglable
  return (any(xbit(Q)[1:k]) || any(zbit(Q)[1:k])) ? :logical : :trivial
end

"Is `v` non-increasing, up to rounding at the level of `tol`?"
nonincreasing(v; tol = 1e-10) = all(v[i+1] <= v[i] + tol for i in 1:length(v)-1)

# -- State type ----------------------------------------------------------------------------
#
# Tests that build a CAMPS state loop over STATE_TYPES. Since Phase 6 (old CliffordMPS code
# deleted, 2026-11-17) it holds only the MagicTensors `QubitCAMPS`; the loops stay so that a
# second implementation can be staged the same way.

import MagicTensors as mt

const STATE_TYPES = (mt.QubitCAMPS,)
stackname(::Type{mt.QubitCAMPS}) = "QubitCAMPS"

"|0…0⟩ with the identity frame."
fresh_state(::Type{mt.QubitCAMPS}, N) = mt.QubitCAMPS(N)

"""
|0…0⟩ behind the frame with stored dagger `Cdag`, i.e. |ψ⟩ = Cdag⁻¹|0…0⟩ (CONVENTIONS.md §4).
`daggered = true` stores C† itself (upstream's lazy flag set), `false` stores C; both are the
same physical state.
"""
function state_with_frame_dagger(::Type{mt.QubitCAMPS}, Cdag; daggered = false)
  N = nqubits(Cdag)
  mps = ITensorMPS.MPS(ComplexF64, ITensorMPS.siteinds("Qubit", N), "0")
  stored = daggered ? copy(Cdag) : inv(Cdag)
  return mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(stored), daggered)
end

nqubits_state(ψ::mt.QubitCAMPS) = mt.nsites(ψ)

"The MPS part of a CAMPS (not a copy)."
mps_part(ψ::mt.QubitCAMPS) = mt.get_mps(ψ)

"The frame C† as a `CliffordOperator` copy (frame-dependent tests only)."
frame_dagger(ψ::mt.QubitCAMPS) = copy(mt.get_clifford(ψ; daggered = true).clifford)

"""
Q = C†PC, the Heisenberg image of P through the frame, sign included. Computed with
QuantumClifford's `C * Stabilizer` (= C P C†, pinned in TestGoldenSet.jl) applied to C†.
"""
heisenberg_image(ψ, P::PauliOperator) = cliffconj(frame_dagger(ψ), P)

"⟨ψ|P|ψ⟩ of a signed Pauli string on the physical state, not normalised (CONVENTIONS.md §5)."
state_expectation(ψ::mt.QubitCAMPS, P::PauliOperator) = mt.expectation(ψ, P)

"""
Does `QubitCAMPS` take its Clifford branch at angle ϕ (frame-dependent: such gates leave k
unchanged)? Every multiple of π/4 within atol = 1e-10, rtol = 0 (user decision 2026-11-12,
plan "Phase 2 decisions"; the spec leaves the defaults to the port).
"""
recognised_clifford(::Type{mt.QubitCAMPS}, ϕ) = abs(ϕ - round(ϕ / (π / 4)) * (π / 4)) <= 1e-10
