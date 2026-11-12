import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using .DisentangleCAMPS

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
      pushfirst!(digits, n % 0x4)
      n ÷= 0x4
    end
    if ndigits > 0
        prepend!(digits, zeros(Int, max(0, ndigits - length(digits))))
    elseif isempty(digits)
        push!(digits, 0)
    end
    return reverse(digits)
end

function QuantumClifford.PauliOperator(pstr::pp.PauliString)
  coeff_to_phase = Dict(1 => 0x0, im => 0x1, -1 => 0x2, -im => 0x3)
  phase = get(coeff_to_phase, pstr.coeff, nothing)
  isnothing(phase) && throw(ArgumentError("QuantumClifford.PauliOperator only supports ±1 and ±i as phases"))
  return pauliop_from_ppterm(pstr.term, pstr.nqubits; phase = phase)
end

"""```pauliop_from_ppterm(term, N; phase = 0x0)```

The QuantumClifford `PauliOperator` of a PauliPropagation term integer on N qubits, with
the given phase byte. Builds the operator with QuantumClifford's own
`PauliOperator(phase, x, z)` constructor, so it does not depend on the
`PauliOperator(::Vector{Int})` method that CliffordMPS adds (needed by
MagicTensorsPauliPropagation.jl after CliffordMPS is removed).

Layout (test/CONVENTIONS.md §3): qubits are 1-based, j = 1…N, and
`term = Σⱼ codeⱼ · 4^(j−1)` with codes I=0, X=1, Y=2, Z=3. So the code of qubit j is
`(term >> 2(j−1)) & 0b11` (what `pp.getpauli(term, j)` returns): the shift counts bits from
0 at the least significant bit, and qubit 1 sits in the lowest two bits. Qubit j of the
result has x = (code ∈ {X, Y}), z = (code ∈ {Y, Z})."""
function pauliop_from_ppterm(term::Integer, N::Integer; phase::UInt8 = 0x0)
  x = Vector{Bool}(undef, N)
  z = Vector{Bool}(undef, N)
  for j in 1:N
    code = Int(pp.getpauli(term, j))
    x[j] = code == 1 || code == 2
    z[j] = code == 2 || code == 3
  end
  return PauliOperator(phase, x, z)
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

# -- Get leftovers from CAMPS evo ---------------------------------------------

leftover_rotgates(s::Integer, rotations, phases) = rotations[s+1:end], -2 .* phases[s+1:end]

# -- Compute expval of pp.PauliSum on CAMPS -----------------------------------
function cmps.expectation(ψ::cmps.CAMPS, op::pp.PauliSum; verbose = false)
  verbose && println("Converting $(length(op))-term sum…")
  op_cmps, conversiontime, _... = @timed cmps.PauliSum(op)
  verbose && println("Done in $conversiontime s.")
  verbose && println("Computing expectation value of $(length(op))-term sum on CAMPS with bond dims $(ψ.mps)…")
  ev, evtime, _... = @timed cmps.expectation(ψ, op_cmps)
  verbose && println("Done in $evtime s.")
  return ev
end