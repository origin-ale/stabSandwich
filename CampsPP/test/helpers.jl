# Shared helpers for the CampsPP tests. Expected values come from test/CONVENTIONS.md and
# from the dense reference (DisentangleCAMPS/test/DenseReference.jl), never from CampsPP.

using CampsPP
using DisentangleCAMPS
using QuantumClifford
using Random
using Test
import CliffordMPS as cmps
import PauliPropagation as pp

isdefined(@__MODULE__, :DenseReference) ||
  include(joinpath(@__DIR__, "..", "..", "DisentangleCAMPS", "test", "DenseReference.jl"))
const DR = DenseReference

const EXACT_ATOL = 1e-10          # CONVENTIONS.md §6, exact regime

# -- Conversions to the dense reference (CONVENTIONS.md §2, §3) --------------------------

"A `pp.PauliRotation` as an N-character string: its symbols on its qubits, I elsewhere."
function gatestring(rot, N)
  letters = fill('I', N)
  for (sym, q) in zip(rot.symbols, rot.qinds)
    letters[q] = string(sym)[1]
  end
  return String(letters)
end

"PP term integer → string, by the bit rule of §3: code_j = bits 2(j−1), 2(j−1)+1; I,X,Y,Z = 0,1,2,3."
termstring(term, N) = String(["IXYZ"[Int((term >> (2 * (j - 1))) & 0x3) + 1] for j in 1:N])

"⟨ψ|O|ψ⟩ of a `pp.PauliSum` on a dense state."
dense_expect(ψ, O::pp.PauliSum) =
  sum(coeff * DR.expect(ψ, termstring(term, O.nqubits)) for (term, coeff) in O.terms)

"""
Dense ⟨O⟩ after 0 gates and after each gate index in `record_at`, for the circuit
(gates, phases) in the exp(iϕP) convention, starting from the basis state with `ones` in |1⟩.
"""
function dense_dynamics(N, ones, gates, phases, O, record_at)
  ψ = DR.basis_state(N, ones)
  evs = [dense_expect(ψ, O)]
  for (i, (g, ϕ)) in enumerate(zip(gates, phases))
    DR.rotate!(ψ, gatestring(g, N), ϕ)
    i in record_at && push!(evs, dense_expect(ψ, O))
  end
  return evs
end

maxabsdiff(a, b) = length(a) == length(b) ? maximum(abs.(Float64.(a) .- Float64.(b))) : Inf

# -- Expected circuit structure (CONVENTIONS.md §2) --------------------------------------

"Bonds of one bricklayer layer, in gate order: (1,2), (3,4), … then (2,3), (4,5), …"
spec_bonds(N) = vcat([[i, i + 1] for i in 1:2:N-1], [[i, i + 1] for i in 2:2:N-1])

"Layer ends of a t-layer bricklayer XXZ circuit: 3(N−1)·l."
spec_layer_ends(N, t) = [3 * (N - 1) * l for l in 1:t]

"Is `v` non-increasing, up to rounding at the level of `tol`?"
nonincreasing(v; tol = 1e-10) = all(v[i+1] <= v[i] + tol for i in 1:length(v)-1)
