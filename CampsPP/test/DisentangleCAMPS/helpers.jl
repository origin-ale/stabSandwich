# Shared helpers for the DisentangleCAMPS tests (included once from runtests.jl).

using QuantumClifford
using CliffordMPS
using Test
import ITensorMPS

"""
The gate `gate` acting on qubits `inds` of an N-qubit register, as a CliffordOperator.
Contract pinned in TestDisentangler.jl ("assumption pins: controlled gates and embedding").
Relies on the vendored `apply!(::CliffordOperator, gate, sites)`; replace it in the port.
"""
function embed(gate, N, inds)
  D = one(CliffordOperator, N)
  apply!(D, gate, inds)
  return D
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

# -- The two CAMPS stacks (plan, Phase 2: staged side by side until T-post ends) ----------
#
# Tests that loop over STATE_TYPES run on the old CliffordMPS `CAMPS` and on the MagicTensors
# `QubitCAMPS`. When the old code is deleted, drop `CAMPS` here, delete the `CAMPS` methods
# below, and turn every `broken_on_old(T)` branch into a plain `@test`.

import MagicTensors as mt

const STATE_TYPES = (CAMPS, mt.QubitCAMPS)
stackname(::Type{CAMPS}) = "old CAMPS"
stackname(::Type{mt.QubitCAMPS}) = "QubitCAMPS"

"""
Is this a hard-spot-6 special-angle case that the old `apply!` gets wrong (CONVENTIONS.md §2)?
Such checks are `@test_broken` on the old type and plain `@test` on `QubitCAMPS`.
"""
broken_on_old(T) = T === CAMPS

"|0…0⟩ with the identity frame."
fresh_state(::Type{CAMPS}, N) = CAMPS(N)
fresh_state(::Type{mt.QubitCAMPS}, N) = mt.QubitCAMPS(N)

"""
|0…0⟩ behind the frame with stored dagger `Cdag`, i.e. |ψ⟩ = Cdag⁻¹|0…0⟩ (CONVENTIONS.md §4).
For `QubitCAMPS`, `daggered = true` stores C† itself (upstream's lazy flag set), `false`
stores C; both are the same physical state.
"""
state_with_frame_dagger(::Type{CAMPS}, Cdag; daggered = false) = CAMPS(copy(Cdag))
function state_with_frame_dagger(::Type{mt.QubitCAMPS}, Cdag; daggered = false)
  N = nqubits(Cdag)
  mps = ITensorMPS.MPS(ComplexF64, ITensorMPS.siteinds("Qubit", N), "0")
  stored = daggered ? copy(Cdag) : inv(Cdag)
  return mt.QubitCAMPS(mps, mt.QubitCliffordUnitary(stored), daggered)
end

nqubits_state(ψ::CAMPS) = length(ψ)
nqubits_state(ψ::mt.QubitCAMPS) = mt.nsites(ψ)

"The MPS part of a CAMPS (not a copy)."
mps_part(ψ::CAMPS) = ψ.mps
mps_part(ψ::mt.QubitCAMPS) = mt.get_mps(ψ)

"The frame C† as a `CliffordOperator` copy (frame-dependent tests only)."
frame_dagger(ψ::CAMPS) = copy(ψ.Cdag)
frame_dagger(ψ::mt.QubitCAMPS) = copy(mt.get_clifford(ψ; daggered = true).clifford)

"""
Q = C†PC, the Heisenberg image of P through the frame, sign included. Computed with
QuantumClifford's `C * Stabilizer` (= C P C†, pinned in TestGoldenSet.jl) applied to C†.
"""
heisenberg_image(ψ, P::PauliOperator) = (frame_dagger(ψ) * Stabilizer([P]))[1]

"⟨ψ|P|ψ⟩ of a signed Pauli string on the physical state, not normalised (CONVENTIONS.md §5)."
state_expectation(ψ::CAMPS, P::PauliOperator) = expectation(ψ, P)
state_expectation(ψ::mt.QubitCAMPS, P::PauliOperator) = mt.expectation(ψ, P)

"""
Does the implementation of type T take its Clifford branch at angle ϕ (frame-dependent:
such gates leave k unchanged)? Old code: |ϕ| ∈ {0, π/4, π/2, 3π/4, π}, CONVENTIONS.md §2.
`QubitCAMPS`: every multiple of π/4 within atol = 1e-10, rtol = 0 (user decision 2026-11-12,
plan "Phase 2 decisions"; the spec leaves the defaults to the port).
"""
recognised_clifford(::Type{CAMPS}, ϕ) =
  any(isapprox(abs(ϕ), c; atol = 1e-9) for c in (0.0, π/4, π/2, 3π/4, π))
recognised_clifford(::Type{mt.QubitCAMPS}, ϕ) = abs(ϕ - round(ϕ / (π / 4)) * (π / 4)) <= 1e-10
