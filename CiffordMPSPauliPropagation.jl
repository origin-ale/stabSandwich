import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

# -- PauliSum pp -> cmps ------------------------------------------------------

function cmps.PauliSum(pstr::pp.PauliString)
  N = pstr.nqubits
  coeff = pstr.coeff
  pstring_int = pstr.term
  pstring_stab = _inttostabilizer(pstring_int, N)
  return cmps.PauliSum([coeff], pstring_stab)
end

function cmps.PauliSum(psum::pp.PauliSum)
  N = psum.nqubits
  coeffs = collect(values(psum.terms))
  pstrings_int = collect(keys(psum.terms))
  psum_stab = _inttostabilizer(pstrings_int, N)
  return cmps.PauliSum(coeffs, psum_stab)
end

_inttopaulioperator(pstr_int::Integer, nqubits::Integer) = 
cmps.PauliOperator([Int(pp.getpauli(pstr_int, i)) for i in 1:nqubits])

_inttostabilizer(pstr_int::Integer, nqubits::Integer) = 
Stabilizer([_inttopaulioperator(pstr_int, nqubits)])

_inttostabilizer(pstr_ints::Vector{<:Integer}, nqubits::Integer) =
vcat(_inttostabilizer.(pstr_ints, nqubits)...)

# -- pp.PauliString -> QuantumClifford.PauliOperator --------------------------

function paulivec(n::Integer, ndigits::Int=0)
    n >= 0 || throw(ArgumentError("n must be non-negative"))
    digits = Int[]
    while n > 0
        pushfirst!(digits, n % 4)
        n ÷= 4
    end
    if ndigits > 0
        prepend!(digits, zeros(Int, max(0, ndigits - length(digits))))
    elseif isempty(digits)
        push!(digits, 0)
    end
    return reverse(digits)
end

function QuantumClifford.PauliOperator(pstr::pp.PauliString) # TODO: implement with Dict
  coeff = pstr.coeff
  if coeff == 1
    phase = 0x0
  elseif coeff == im
    phase = 0x1
  elseif coeff == -1
    phase = 0x2
  elseif coeff == -im
    phase = 0x3
  else 
    throw(ArgumentError("QuantumClifford.PauliOperator only supports ±1 and ±i as phases"))
  end
  return PauliOperator(paulivec(pstr.term, pstr.nqubits); phase = phase)
end

# -- Extract pp.PauliStrings from pp.PauliRotations ---------------------------

getpauli(rot::pp.PauliRotation, N::Integer) = pp.PauliString(N, rot.symbols, rot.qinds)

# -- Evolving on pp objects ---------------------------------------------------

DisentangleCAMPS.evolve_bonddim(ψ::cmps.CAMPS, 
                        χ::Integer, 
                        paulistrings::Vector{<:pp.PauliString},
                        phases::Vector{<:Real}; 
                        showprogress = false) = 
DisentangleCAMPS.evolve_bonddim(ψ, χ, PauliOperator.(paulistrings), phases; showprogress = showprogress)

DisentangleCAMPS.evolve_bonddim(ψ::cmps.CAMPS, 
                        χ::Integer, 
                        paulirots::Vector{<:pp.PauliRotation},
                        phases::Vector{<:Real}; 
                        showprogress = false) = 
DisentangleCAMPS.evolve_bonddim(ψ, χ, getpauli.(paulirots, length(ψ)), phases; showprogress = showprogress)

DisentangleCAMPS.evolve(ψ::cmps.CAMPS, 
                        t::Integer, 
                        paulistrings::Vector{<:pp.PauliString},
                        phases::Vector{<:Real}; 
                        showprogress = false) = 
evolve(ψ, t, PauliOperator.(paulistrings), phases; showprogress = showprogress)

DisentangleCAMPS.evolve(ψ::cmps.CAMPS, 
                        t::Integer, 
                        paulirots::Vector{<:pp.PauliRotation},
                        phases::Vector{<:Real}; 
                        showprogress = false) = 
evolve(ψ, t, getpauli.(paulirots, length(ψ)), phases; showprogress = showprogress)