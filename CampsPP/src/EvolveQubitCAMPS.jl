# MagicTensors `QubitCAMPS` versions of the functions in Evolve.jl.
#
# Staged side by side with the old CliffordMPS `CAMPS` methods (plan, Phase 2) and selected by
# dispatch on the state type; the old methods stay unchanged until T-post, then this file is
# the only implementation.
#
# Frame convention (test/CONVENTIONS.md §4): `QubitCAMPS` stores the state-frame Clifford C,
# |ψ⟩ = C |mps⟩, with a lazy `daggered` flag. Every frame write below goes through C†, fetched
# with `get_clifford(ψ; daggered = true)` immediately before the write: a reference held across
# a lazy flip is stale, and writes to it are silently lost.
#
# Never call QuantumClifford's 2-argument `apply!(::CliffordOperator, op)` here: it defaults to
# `phases = false`. Products of `CliffordOperator`s use QuantumClifford's native `A * B`, the
# operator product (B acts first).

import MagicTensors as mt
using Random: AbstractRNG

"Default absolute tolerance for recognising a Clifford angle (a multiple of π/4)."
const CLIFFORD_ATOL = 1e-10
"Default relative tolerance for recognising a Clifford angle (a multiple of π/4)."
const CLIFFORD_RTOL = 0.0

bonddim(ψ::mt.QubitCAMPS) = maximum(dim.(linkinds(mt.get_mps(ψ))); init = 1)

# == Gates ==============================================================================

"""```rotation_gate(P, ϕ)```

The rotation exp(iϕP) = cos ϕ · 1 + i sin ϕ · P as a two-term `MagicTensors.QubitPauliSum`
(test/CONVENTIONS.md §2). The sign of the Hermitian string `P` is folded into its coefficient.
`P` is copied, because a `QubitPauliSum` keeps the key vector it is given."""
function rotation_gate(P::PauliOperator, ϕ::Real)
  check_hermitian(P)
  N = nqubits(P)
  R = mt.QubitPauliSum(N)
  push!(R, mt.QubitPauli(N), cos(ϕ))
  push!(R, mt.QubitPauli(copy(P)), im * sin(ϕ))
  return R
end

"A gate string must be Hermitian (phase byte 0 or 2, i.e. sign ±1); ±i would corrupt the frame."
check_hermitian(P::PauliOperator) = iseven(P.phase[]) ||
  throw(ArgumentError("rotation needs a Hermitian Pauli string (sign ±1), got $P"))

"""
If ϕ is a multiple mπ/4 within `isapprox(ϕ, mπ/4; atol, rtol)`, return m mod 4, the number of
quarter turns of exp(iϕP) up to a global phase; otherwise return `nothing`. Every multiple of
π/4 is recognised as Clifford, including 0, ±π/2 and ±π, so all of them take the exact
frame-update path.
"""
function clifford_quarter_turns(ϕ::Real; atol::Real = CLIFFORD_ATOL, rtol::Real = CLIFFORD_RTOL)
  isfinite(ϕ) || throw(ArgumentError("rotation angle must be finite, got $ϕ"))
  m = round(Int, ϕ / (π / 4))
  isapprox(ϕ, m * π / 4; atol = atol, rtol = rtol) || return nothing
  return mod(m, 4)
end

"""
Tableau of exp(i·n·π/4·P), up to a global phase, for n ∈ 0:3. A Pauli A that anticommutes with
P maps to iPA (n = 1), −A (n = 2) or −iPA (n = 3); a commuting A is fixed. The sign of P counts.
"""
function pauli_rotation_tableau(P::PauliOperator, n::Integer)
  N = nqubits(P)
  G = one(CliffordOperator, N)
  for i in 1:2N
    A = G[i]
    if comm(P, A) == 0x1
      G[i] = n == 1 ? (im * P * A) : n == 2 ? -A : ((-im) * P * A)
    end
  end
  return G
end

"Frame update for the Clifford rotation exp(i·n·π/4·P): C ← G·C, i.e. C† ← C†·G†."
function apply_clifford_rotation!(ψ::mt.QubitCAMPS, P::PauliOperator, n::Integer)
  n == 0 && return nothing
  Gdag = pauli_rotation_tableau(P, mod(-n, 4))
  Cdag = mt.get_clifford(ψ; daggered = true)
  Cdag.clifford = Cdag.clifford * Gdag
  return nothing
end

"The Heisenberg image C†PC of `P` through the frame of `ψ`, sign included."
heisenberg_image(ψ::mt.QubitCAMPS, P::PauliOperator) =
  mt.conjugate(mt.QubitPauli(P), mt.get_clifford(ψ; daggered = true)).pauli

# == apply! =============================================================================

"""```apply!(ψ::QubitCAMPS, k, P, ϕ; atol = $CLIFFORD_ATOL, rtol = $CLIFFORD_RTOL)```

Apply exp(iϕP) to the CAMPS ψ with k magic qubits, disentangling if possible.
Modify ψ in place and return the new number of magic qubits.

An angle within `isapprox(ϕ, mπ/4; atol, rtol)` of a multiple of π/4 is applied exactly to the
Clifford frame and leaves k unchanged. Otherwise, with Q = C†PC (`paulinature`):
- `:disentanglable`: C ← C·D and exp(i·sϕ·X_{k+1}) on the MPS, with (D, s) from `disentangler`;
  k grows by one;
- `:logical`: exp(iϕQ) on the MPS (MPO apply, ITensorMPS default cutoff);
- `:trivial`: nothing."""
function QuantumClifford.apply!(ψ::mt.QubitCAMPS,
                                k::Integer,
                                P::PauliOperator,
                                ϕ::Real;
                                atol::Real = CLIFFORD_ATOL,
                                rtol::Real = CLIFFORD_RTOL)
  check_hermitian(P)
  n = clifford_quarter_turns(ϕ; atol = atol, rtol = rtol)
  if !isnothing(n)
    apply_clifford_rotation!(ψ, P, n)
    return k
  end

  Q = heisenberg_image(ψ, P)
  nature = paulinature(k, Q)
  if nature == :disentanglable
    D, sign = disentangler(k, Q)
    mt.apply_to_clifford_dagger!(ψ, mt.QubitCliffordUnitary(inv(D))) # C ← C·D
    addmagicstate!(ψ, k, sign * ϕ)
    k += 1
  elseif nature == :logical
    mt.apply_to_mps!(ψ, rotation_gate(Q, ϕ))
  end
  return k
end

"""```apply!(ψ::QubitCAMPS, P, ϕ, thl; atol = $CLIFFORD_ATOL, rtol = $CLIFFORD_RTOL)```

Apply exp(iϕP) to the CAMPS ψ, keeping the exact Clifford handling of `apply!(ψ, k, P, ϕ)` but
folding every non-Clifford gate directly into the MPS part, exp(iϕ·C†PC), with SVD truncation
at threshold `thl` (no disentangling, no new magic qubits). Modify ψ in place, return nothing."""
function QuantumClifford.apply!(ψ::mt.QubitCAMPS,
                                P::PauliOperator,
                                ϕ::Real,
                                thl::Real;
                                atol::Real = CLIFFORD_ATOL,
                                rtol::Real = CLIFFORD_RTOL)
  check_hermitian(P)
  n = clifford_quarter_turns(ϕ; atol = atol, rtol = rtol)
  if !isnothing(n)
    apply_clifford_rotation!(ψ, P, n)
    return nothing
  end
  mt.apply_to_mps!(ψ, rotation_gate(heisenberg_image(ψ, P), ϕ); cutoff = thl)
  return nothing
end

"""```addmagicstate!(ψ::QubitCAMPS, k, α)```

Multiply qubit k+1 of the MPS of ψ by exp(iαX) = cos α · 1 + i sin α · X, turning its |0⟩ into
the Liu and Clark (2025) magic state. Exact single-site multiply: unlike the old MPO route it
does not re-gauge or compress the MPS (same state; see the plan's parity notes)."""
function addmagicstate!(ψ::mt.QubitCAMPS, k::Integer, α::Real)
  c, s = cos(α), sin(α)
  mt.apply_to_mps!(ψ, ComplexF64[c im*s; im*s c], [k + 1])
  return nothing
end

# == Circuits ===========================================================================

"""```evolve_bonddim(ψ::QubitCAMPS, χ, paulis, phases; [showprogress, atol, rtol])```

Evolve the CAMPS ψ along the Pauli rotation circuit with Pauli strings `paulis` and the given
phases, until the end or bond dimension χ. Return the evolved CAMPS, k and the stopping gate."""
function evolve_bonddim(ψ::mt.QubitCAMPS,
                        χ::Integer,
                        paulistrings::Vector{<:PauliOperator},
                        phases::Vector{<:Real};
                        showprogress = false,
                        atol::Real = CLIFFORD_ATOL,
                        rtol::Real = CLIFFORD_RTOL)
  k = 0
  s = 0
  progressthresh = ProgressUnknown(0; dt = 0.05, desc = "Evolving CAMPS… t =", enabled = showprogress)
  while s < length(paulistrings) && bonddim(ψ) < χ
    s += 1
    k = apply!(ψ, k, paulistrings[s], phases[s]; atol = atol, rtol = rtol)
    next!(progressthresh; showvalues = generate_showvalues(χ, bonddim(ψ)))
  end
  finish!(progressthresh)
  return ψ, k, s
end

"""```evolve(ψ::QubitCAMPS, t, paulis, phases; [showprogress, atol, rtol])```

Evolve the CAMPS ψ along the first t gates of the Pauli rotation circuit with Pauli strings
`paulis` and the given phases. Return the evolved CAMPS and k."""
function evolve(ψ::mt.QubitCAMPS,
                t::Integer,
                paulistrings::Vector{<:PauliOperator},
                phases::Vector{<:Real};
                showprogress = false,
                atol::Real = CLIFFORD_ATOL,
                rtol::Real = CLIFFORD_RTOL)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    k = apply!(ψ, k, paulistrings[s], phases[s]; atol = atol, rtol = rtol)
    next!(progressbar)
  end
  return ψ, k
end

"""```evolve_deepcliffords(rng, ψ::QubitCAMPS, t, paulis, phases; [showprogress, atol, rtol])```

Evolve the CAMPS ψ along a Pauli rotation-doped deep Clifford circuit, until layer t.
Layer s applies a uniformly random N-qubit Clifford R (`QuantumClifford.random_clifford(rng, N)`,
replacing the old 2N²-deep random circuit), then the rotation (paulis[s], phases[s]).
Return the evolved CAMPS and k."""
function evolve_deepcliffords(rng::AbstractRNG,
                              ψ::mt.QubitCAMPS,
                              t::Integer,
                              paulistrings::Vector{<:PauliOperator},
                              phases::Vector{<:Real};
                              showprogress = false,
                              atol::Real = CLIFFORD_ATOL,
                              rtol::Real = CLIFFORD_RTOL)
  N = mt.nsites(ψ)
  k = 0
  progressbar = Progress(t; desc = "Evolving…", enabled = showprogress)
  for s in 1:t
    R = random_clifford(rng, N)
    Cdag = mt.get_clifford(ψ; daggered = true)
    Cdag.clifford = Cdag.clifford * inv(R) # C ← R·C
    k = apply!(ψ, k, paulistrings[s], phases[s]; atol = atol, rtol = rtol)
    next!(progressbar)
  end
  return ψ, k
end
