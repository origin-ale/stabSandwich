# Shared helpers for the DisentangleCAMPS tests (included once from runtests.jl).

using QuantumClifford
using CliffordMPS
using Test

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
