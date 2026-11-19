# Bridge between PauliPropagation and MagicTensors: the pp → `QubitPauliSum` converter, the
# `expectation` seam for `pp.PauliSum` observables on a `QubitCAMPS`, and the evolve wrappers
# taking pp gates, plus the pp → QuantumClifford helpers they use (formerly in
# CliffordMPSPauliPropagation.jl, deleted with the old CliffordMPS stack in Phase 6).

import MagicTensors as mt
import PauliPropagation as pp
import ITensors, ITensorMPS
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

"""```expectation(ψ::QubitCAMPS, op::pp.PauliSum; method = :perterm, normalized = false, verbose = false)```

⟨ψ|O|ψ⟩ = Σ c_P ⟨mps|C†PC|mps⟩ of a PauliPropagation sum on the physical state
|ψ⟩ = C|mps⟩, as a `ComplexF64` (test/CONVENTIONS.md §5, "The pp seam"). This is the single
seam through which every pp observable is evaluated on a `QubitCAMPS`.

- `method = :perterm` (default): convert, conjugate into the frame, then evaluate the terms
  one by one on a copy of the MPS (`perterm_expectation`).
- `method = :mpo`: convert, conjugate, then MagicTensors' OpSum → MPO → `inner`; kept for checks.
- `normalized = true` divides by ⟨mps|mps⟩. An empty sum gives 0.

ψ is left bitwise unchanged (MPS tensors, stored frame and its `daggered` flag): upstream's
`get_clifford(ψ; daggered = true)` would invert a stored C in place, so C† is built here
without touching ψ (`frame_dagger`). Serial: no threads inside one call."""
function mt.expectation(ψ::mt.QubitCAMPS, op::pp.PauliSum;
                        method = :perterm, normalized::Bool = false, verbose = false)
  method in (:perterm, :mpo) ||
    throw(ArgumentError("unknown method $(repr(method)); use :perterm or :mpo"))
  verbose && println("Converting $(length(op))-term sum…")
  op_mt, conversiontime, _... = @timed mt.QubitPauliSum(op)
  verbose && println("Done in $conversiontime s.")
  mps = mt.get_mps(ψ)
  verbose && println("Computing expectation value of $(length(op))-term sum ($method) on CAMPS with bond dims $(ITensorMPS.linkdims(mps))…")
  ev, evtime, _... = @timed begin
    op_conj = mt.conjugate(op_mt, frame_dagger(ψ))     # C† P C, signs folded into the coefficients
    if method === :perterm
      perterm_expectation(mps, op_conj; normalized)
    elseif length(op_conj) == 0                          # upstream's MPO route throws on an empty sum
      zero(ComplexF64)
    else
      ComplexF64(mt.expectation(mps, op_conj; normalized))
    end
  end
  verbose && println("Done in $evtime s.")
  return ev
end

"""```frame_dagger(ψ::QubitCAMPS)```

The frame C† of ψ as a `QubitCliffordUnitary`, without changing ψ: the stored unitary itself
when ψ stores C† (read only), else a fresh inverse of the stored C."""
frame_dagger(ψ::mt.QubitCAMPS) = ψ.daggered ? ψ.clifford : inv(ψ.clifford)

const PERTERM_OPNAMES = ("Id", "X", "Y", "Z")   # indexed by code + 1, codes I=0, X=1, Y=2, Z=3

"""```perterm_expectation(mps, Q::QubitPauliSum; normalized = false, nsteps = nothing)```

Σ c_Q ⟨mps|Q|mps⟩ over the (already frame-conjugated) terms of `Q`, term by term, on a copy of
`mps`; `mps` itself is not changed. The copy is put in canonical form with centre at site 1,
so a term whose last non-identity site is m needs only sites 1…m: the right part closes as
the identity on bond m. Terms are sorted by their Pauli codes, and the left environments are
reused along each term's common prefix with the previous one. With `normalized = true` the
sum is divided by ⟨mps|mps⟩. An empty sum gives 0. `Q` must be on `length(mps)` qubits
(`DimensionMismatch` otherwise).

Cost hook (test/CONVENTIONS.md §5, "Cost of `method = :perterm`"): with `nsteps::Ref{Int}`,
the number of left environments built (one per site absorbed into an environment) is added
to `nsteps[]`; the returned value is the same as without it."""
function perterm_expectation(mps::ITensorMPS.MPS, Q::mt.QubitPauliSum; normalized::Bool = false,
                             nsteps::Union{Nothing, Ref{Int}} = nothing)
  N = length(mps)
  mt.nsites(Q) == N ||
    throw(DimensionMismatch("sum on $(mt.nsites(Q)) qubits, MPS on $N sites"))
  terms = Vector{Tuple{Vector{Int8}, ComplexF64}}(undef, length(Q.paulis))
  for (i, (xz, c)) in enumerate(Q.paulis)
    P = PauliOperator(0x00, N, xz)                       # sign already in c (MagicTensors push!)
    terms[i] = ([Int8(xzcode(P[j]...)) for j in 1:N], ComplexF64(c))
  end

  m1 = ITensorMPS.orthogonalize!(deepcopy(mps), 1)
  sites = ITensorMPS.siteinds(m1)
  ket = [m1[j] for j in 1:N]
  bra = [ITensors.dag(ITensors.prime(m1[j])) for j in 1:N]
  ops = [Tuple(ITensors.op(n, sites[j]) for n in PERTERM_OPNAMES) for j in 1:N]
  links = [ITensorMPS.linkind(m1, b) for b in 1:N-1]
  normsq = ComplexF64(ITensors.scalar(ket[1] * ITensors.dag(ket[1])))

  # Close a left environment L over sites 1…m: the trace over bond m (identity on the right).
  close_env(L, m) = m == 0 ? normsq :
                    m == N ? ComplexF64(ITensors.scalar(L)) :
                    ComplexF64(ITensors.scalar(L * ITensors.delta(links[m], ITensors.prime(links[m]))))

  envs = Vector{ITensors.ITensor}(undef, N)              # envs[j]: environment over sites 1…j of prev
  prev = fill(Int8(-1), N)
  valid = 0                                              # envs[1:valid] belong to prev[1:valid]
  ev = zero(ComplexF64)
  for i in sortperm(terms; by = first)
    codes, c = terms[i]
    m = something(findlast(!iszero, codes), 0)           # last non-identity site
    j0 = 0
    while j0 < min(valid, m) && codes[j0 + 1] == prev[j0 + 1]
      j0 += 1
    end
    L = j0 == 0 ? ITensors.ITensor(1.0) : envs[j0]
    for j in j0+1:m
      L = ((L * ket[j]) * ops[j][codes[j] + 1]) * bra[j]
      envs[j] = L
      prev[j] = codes[j]
      isnothing(nsteps) || (nsteps[] += 1)               # counts the contractions done
    end
    valid = m                                            # entries past m belong to an older prefix
    ev += c * close_env(L, m)
  end
  return normalized ? ev / normsq : ev
end

"Code of a qubit's Pauli from its (x, z) bits: I=0, X=1, Y=2, Z=3 (test/CONVENTIONS.md §3)."
xzcode(x::Bool, z::Bool) = x ? (z ? 2 : 1) : (z ? 3 : 0)

"⟨ψ|P|ψ⟩ of a signed QuantumClifford Pauli string on the physical state, not normalised."
mt.expectation(ψ::mt.QubitCAMPS, P::PauliOperator) = mt.expectation(ψ, mt.QubitPauli(P))

# -- Evolving on pp objects ------------------------------------------------------------

# `evolve_bonddim` and `evolve` also accept the circuit as pp Pauli strings or pp Pauli
# rotations (one per gate); they are converted to signed QuantumClifford strings first.

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

# -- pp.PauliString -> QuantumClifford.PauliOperator ----------------------------------

"""```PauliOperator(pstr::pp.PauliString)```

The QuantumClifford `PauliOperator` of a PauliPropagation string, with its coefficient as \
the phase; the coefficient must be exactly ±1 or ±i."""
function QuantumClifford.PauliOperator(pstr::pp.PauliString)
  coeff_to_phase = Dict(1 => 0x0, im => 0x1, -1 => 0x2, -im => 0x3)
  phase = get(coeff_to_phase, pstr.coeff, nothing)
  isnothing(phase) && throw(ArgumentError("QuantumClifford.PauliOperator only supports ±1 and ±i as phases"))
  return pauliop_from_ppterm(pstr.term, pstr.nqubits; phase = phase)
end

"""```pauliop_from_ppterm(term, N; phase = 0x0)```

The QuantumClifford `PauliOperator` of a PauliPropagation term integer on N qubits, with
the given phase byte, built with QuantumClifford's own `PauliOperator(phase, x, z)`
constructor.

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

# -- Extract pp.PauliStrings from pp.PauliRotations ------------------------------------

"""```getpauli(rot::pp.PauliRotation, N)```

The N-qubit `pp.PauliString` (coefficient 1) of the generator of the rotation ```rot```. \
Wrap it in ```PauliOperator``` for the CAMPS ```apply!```."""
getpauli(rot::pp.PauliRotation, N::Integer) = pp.PauliString(N, rot.symbols, rot.qinds)

# -- Leftover gates after a CAMPS evolution up to gate s ---------------------------------

"""```leftover_rotgates(s, rotations, phases)```

The gates after gate ```s``` and their PauliPropagation angles −2ϕ (```rotations``` and \
```phases``` as built by ```rotation_circuit``` or ```xxz_circuit```, gate exp(iϕP)), ready \
for ```pp.propagate```."""
leftover_rotgates(s::Integer, rotations, phases) = rotations[s+1:end], -2 .* phases[s+1:end]
