# Bridge between PauliPropagation and MagicTensors: the pp → `QubitPauliSum` converter, the
# `expectation` seam for `pp.PauliSum` observables on a `QubitCAMPS`, and the evolve wrappers
# taking pp gates. Staged next to CliffordMPSPauliPropagation.jl (plan, Phase 2), which keeps
# the old CliffordMPS versions until T-post.

import MagicTensors as mt
import PauliPropagation as pp
using QuantumClifford
using .DisentangleCAMPS

# -- pp.PauliSum -> MagicTensors.QubitPauliSum -----------------------------------------

"""```QubitPauliSum(psum::pp.PauliSum)```

Convert a PauliPropagation sum to a `MagicTensors.QubitPauliSum` term by term: same strings
(bit rule of test/CONVENTIONS.md §3), coefficients converted to `ComplexF64`. No splatting,
so the cost is linear in the number of terms and nothing is compiled per sum size."""
function mt.QubitPauliSum(psum::pp.PauliSum)
  N = psum.nqubits
  out = mt.QubitPauliSum(N)
  for (term, coeff) in psum.terms
    # A fresh key vector per term: the sum keeps the vector it is given.
    push!(out, mt.QubitPauli(pauliop_from_ppterm(term, N)), ComplexF64(coeff))
  end
  return out
end

"""```QubitPauliSum(pstr::pp.PauliString)```

One-term `MagicTensors.QubitPauliSum` with the string's coefficient as `ComplexF64`."""
function mt.QubitPauliSum(pstr::pp.PauliString)
  out = mt.QubitPauliSum(pstr.nqubits)
  push!(out, mt.QubitPauli(pauliop_from_ppterm(pstr.term, pstr.nqubits)), ComplexF64(pstr.coeff))
  return out
end

# -- Expectation values on QubitCAMPS --------------------------------------------------

"""```expectation(ψ::QubitCAMPS, op::pp.PauliSum; verbose = false)```

⟨ψ|O|ψ⟩ of a PauliPropagation sum on the physical state |ψ⟩ = C|mps⟩, not normalised
(test/CONVENTIONS.md §5). This is the single seam through which every pp observable is
evaluated on a `QubitCAMPS`: convert, then MagicTensors' conjugate → OpSum → MPO → inner."""
function mt.expectation(ψ::mt.QubitCAMPS, op::pp.PauliSum; verbose = false)
  verbose && println("Converting $(length(op))-term sum…")
  op_mt, conversiontime, _... = @timed mt.QubitPauliSum(op)
  verbose && println("Done in $conversiontime s.")
  verbose && println("Computing expectation value of $(length(op))-term sum on CAMPS with bond dims $(linkdims(mt.get_mps(ψ)))…")
  ev, evtime, _... = @timed mt.expectation(ψ, op_mt)
  verbose && println("Done in $evtime s.")
  return ev
end

"⟨ψ|P|ψ⟩ of a signed QuantumClifford Pauli string on the physical state, not normalised."
mt.expectation(ψ::mt.QubitCAMPS, P::PauliOperator) = mt.expectation(ψ, mt.QubitPauli(P))

# -- Evolving on pp objects ------------------------------------------------------------

DisentangleCAMPS.evolve_bonddim(ψ::mt.QubitCAMPS,
                                χ::Integer,
                                paulistrings::Vector{<:pp.PauliString},
                                phases::Vector{<:Real};
                                kwargs...) =
DisentangleCAMPS.evolve_bonddim(ψ, χ, PauliOperator.(paulistrings), phases; kwargs...)

DisentangleCAMPS.evolve_bonddim(ψ::mt.QubitCAMPS,
                                χ::Integer,
                                paulirots::Vector{<:pp.PauliRotation},
                                phases::Vector{<:Real};
                                kwargs...) =
DisentangleCAMPS.evolve_bonddim(ψ, χ, PauliOperator.(getpauli.(paulirots, mt.nsites(ψ))), phases; kwargs...)

DisentangleCAMPS.evolve(ψ::mt.QubitCAMPS,
                        t::Integer,
                        paulistrings::Vector{<:pp.PauliString},
                        phases::Vector{<:Real};
                        kwargs...) =
DisentangleCAMPS.evolve(ψ, t, PauliOperator.(paulistrings), phases; kwargs...)

DisentangleCAMPS.evolve(ψ::mt.QubitCAMPS,
                        t::Integer,
                        paulirots::Vector{<:pp.PauliRotation},
                        phases::Vector{<:Real};
                        kwargs...) =
DisentangleCAMPS.evolve(ψ, t, PauliOperator.(getpauli.(paulirots, mt.nsites(ψ))), phases; kwargs...)
