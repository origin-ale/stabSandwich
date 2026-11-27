# Shared helpers for the CampsPP tests. Expected values come from test/CONVENTIONS.md and
# from the dense reference (CampsPP/test/DisentangleCAMPS/DenseReference.jl), never from CampsPP.

using CampsPP
using CampsPP.DisentangleCAMPS
using QuantumClifford
using Random
using Test
import MagicTensors as mt
import PauliPropagation as pp

isdefined(@__MODULE__, :DenseReference) ||
  include(joinpath(@__DIR__, "DisentangleCAMPS", "DenseReference.jl"))
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

# Kicked Ising (CONVENTIONS.md §2 "Kicked Ising"): the circuit and its dense profile from the
# spec, in KimSpec.jl. (The KIM tests were written before the code, KIM plan step 2, and their
# gates were removed in step 4.)
include(joinpath(@__DIR__, "KimSpec.jl"))

"Is `v` non-increasing, up to rounding at the level of `tol`?"
nonincreasing(v; tol = 1e-10) = all(v[i+1] <= v[i] + tol for i in 1:length(v)-1)

# -- State types -------------------------------------------------------------------------

# Tests that build a CAMPS state loop over STATE_TYPES. Since Phase 6 (old CliffordMPS code
# deleted) it holds only the MagicTensors `QubitCAMPS`; the loops stay so that a second
# implementation can be staged the same way.
const STATE_TYPES = (mt.QubitCAMPS,)
stackname(::Type{mt.QubitCAMPS}) = "QubitCAMPS"

"⟨ψ|O|ψ⟩ of a `pp.PauliSum` on the physical state, not normalised (CONVENTIONS.md §5)."
state_expectation(ψ::mt.QubitCAMPS, O::pp.PauliSum) = mt.expectation(ψ, O)
"⟨ψ|P|ψ⟩ of a signed QuantumClifford Pauli string on the physical state, not normalised."
state_expectation(ψ::mt.QubitCAMPS, P::PauliOperator) = mt.expectation(ψ, P)

nqubits_state(ψ::mt.QubitCAMPS) = mt.nsites(ψ)

# Cost hook (CONVENTIONS.md §5, "Cost of `method = :perterm`"), used by TestPerTermSeam.jl:
# `CampsPP.perterm_expectation(mps, Q::QubitPauliSum; nsteps::Ref{Int})` adds to `nsteps[]` the
# number of left environments it builds, one per site absorbed into an environment, and returns
# the same value as without `nsteps`. A sum whose number of qubits differs from `length(mps)`
# raises `DimensionMismatch` (review fix, 2026-11-19). The item 4 tests used to be gated on this
# hook and on the seam's `method`/`normalized` keywords; the gates were removed once item 4 landed.

"The frame C† as a `CliffordOperator` copy (frame-dependent tests only)."
frame_dagger(ψ::mt.QubitCAMPS) = copy(mt.get_clifford(ψ; daggered = true).clifford)

# -- Child processes ----------------------------------------------------------------------

const CHILD_HEADER = raw"""
report(name, value) = println("RESULT\t", name, "\t", repr(value))
function errorinfo(f)
  try
    f()
    return (false, "returned normally")
  catch e
    return (e isa ErrorException, sprint(showerror, e))
  end
end
datarows(path) = [parse(Int, split(l, '\t')[1]) for l in eachline(path)
                  if !isempty(strip(l)) && !startswith(l, '#')]
"""

"""
Run `code` as a script in a fresh Julia process (this binary and its flags, then `flags`, no
startup file, the active test environment) and return the values it reported with
`report(name, value)`, as a Dict name => value. A child that fails makes this throw.
"""
function run_child(code::AbstractString; flags::Cmd = ``)
  out = mktempdir() do dir
    script = joinpath(dir, "child.jl")
    write(script, CHILD_HEADER * code)
    read(`$(Base.julia_cmd()) $flags --startup-file=no --project=$(Base.active_project()) $script`, String)
  end
  results = Dict{String, Any}()
  for line in eachline(IOBuffer(out))
    startswith(line, "RESULT\t") || continue
    _, name, value = split(line, '\t'; limit = 3)
    results[name] = Core.eval(@__MODULE__, Meta.parse(value))
  end
  return results
end
