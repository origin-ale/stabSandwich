# Layer-end disentangler for MagicTensors `QubitCAMPS`: a local reimplementation of the old
# CliffordMPS `disentangle!(camps, strategy, max_iter; criterion, min_diff, svd_kwargs...)`
# (plan, Phase 3), built on the upstream disentangler interface. Upstream has neither the
# `:chi3` criterion, nor `min_diff` stopping, nor the old gate set and sweep orders.
#
# What is reproduced from the old code, so that both stacks make the same choices:
# - the gate set: the old `Clifford2IndexSet(:entangle)`, 19 two-qubit Cliffords in the old
#   order (upstream's `QubitCliffordGateSet(2, :entangle)` is a different set);
# - the bond search: every gate is tried on a copy of the MPS orthogonalized at the bond, the
#   first minimum of the criterion wins, and it is applied only if strictly below the identity;
#   the real MPS is touched only when a gate is applied;
# - the criteria: `:chi3` = (number of singular values)³ and `:entangle` = von Neumann entropy
#   in bits of the squared singular values;
# - the sweeps: `:snake` (forward then backward over the bonds, skipping bonds of dimension 1 in
#   a snapshot taken at the start of each direction) and `:brickwork` (odd bonds, then even
#   bonds in reverse);
# - the stopping rule: repeat the sweep up to `max_iter` times, stop once its total gain is
#   ≤ `min_diff`.
#
# `svd_kwargs` go only to the SVD of the applied gate; the search SVDs are untruncated, as in
# the old code (its `disentangler(mps, nL, nR)` does not pass them on for neighbouring sites).
# ITensors keeps exact zero singular values, so in the search the number of singular values is
# min(2χ_L, 2χ_R) for every candidate gate: on these nearest-neighbour sweeps `:chi3` never
# finds an improvement and never applies a gate, in the old code and here alike.
#
# Frame writes go through `mt.transform!`, which fetches C† immediately before writing it.

import MagicTensors as mt

# == Gate set ===========================================================================

# Old one-qubit factors stA[1 + a] for a = 0:3 (only these occur in the set): I, H, S, H·S.
const _OLD_A = (tId1, tHadamard, tPhase, tHadamard * tPhase)
const _OLD_ISWAP = (tId1 ⊗ tHadamard) * tSWAP * tCNOT * tSWAP * tCNOT * (tHadamard ⊗ tId1) *
                   (tPhase ⊗ tPhase)
# Old entangling cores stR[1 + r] for r = 1:3: SWAP, iSWAP, CNOT.
const _OLD_R = (tSWAP, _OLD_ISWAP, tCNOT)

# (r, a1, a2) of the old `entangle_set` numbers 576, 1152, …, 9000, in order. All have
# b1 = b2 = 0 and s1 = s2 = 0, so the old gate is stR[1 + r] * (stA[1 + a1] ⊗ stA[1 + a2]).
const _OLD_ENTANGLE_SET = ((1, 0, 0),
  (2, 0, 0), (2, 0, 1), (2, 0, 3), (2, 1, 0), (2, 1, 1), (2, 1, 3), (2, 3, 0), (2, 3, 1), (2, 3, 3),
  (3, 0, 0), (3, 0, 1), (3, 0, 2), (3, 1, 0), (3, 1, 1), (3, 1, 2), (3, 3, 0), (3, 3, 1), (3, 3, 2))

"""```old_entangle_gates()```

The 19 two-qubit Cliffords of the old CliffordMPS `Clifford2IndexSet(:entangle)`, in the old
order, as `QuantumClifford.CliffordOperator`s."""
old_entangle_gates() =
  [_OLD_R[r] * (_OLD_A[1 + a1] ⊗ _OLD_A[1 + a2]) for (r, a1, a2) in _OLD_ENTANGLE_SET]

# == Criteria ===========================================================================

# Score of the singular values `spec` across the bond; lower is less entangled.
criterion_value(::Val{:chi3}, spec) = float(length(spec)^3)
function criterion_value(::Val{:entangle}, spec)
  ee = 0.0 # Same summation order as the old `eEntropy`
  for σ in spec
    p = σ^2
    if p > 0
      ee -= p * log2(p)
    end
  end
  return ee
end

# == Bond search ========================================================================

"""```BondSearchDisentangler(criterion; svd_kwargs...)```

`MagicTensors.AbstractLocalDisentangler` with the old CliffordMPS bond search: on bond
(n, n+1), try each gate of `old_entangle_gates()`, score the resulting singular values with
`criterion` (`:chi3` or `:entangle`), and apply the first best gate if it scores strictly below
the identity. The search SVDs are untruncated; `svd_kwargs` apply to the applied gate only.
`disentangle!(de, ψ, n, move_ortho_right)` returns `(1, gain, gate index)` or `(0, 0.0, 0)`."""
struct BondSearchDisentangler{C} <: mt.AbstractLocalDisentangler
  gates::Vector{mt.QubitCliffordGate}
  matrices::Vector{Matrix{ComplexF64}}
  svd_kwargs::Dict{Symbol, Any}
end

function BondSearchDisentangler(criterion::Symbol; svd_kwargs...)
  criterion in (:chi3, :entangle) ||
    throw(ArgumentError("criterion must be :chi3 or :entangle, got :$criterion"))
  gates = [mt.QubitCliffordGate(G) for G in old_entangle_gates()]
  return BondSearchDisentangler{criterion}(
    gates, [Matrix{ComplexF64}(Matrix(g)) for g in gates], Dict{Symbol, Any}(svd_kwargs))
end

"Two-qubit identity, scored first as the reference (the old `Clifford2Index(0)`)."
const _ID_2Q = ComplexF64[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]

function bond_value(::BondSearchDisentangler{C}, mps, matrix, n) where {C}
  _, S, _ = mt.apply_and_svd(mps, matrix, n) # Untruncated, as the old search
  return criterion_value(Val(C), diag(S))
end

function mt.disentangle!(de::BondSearchDisentangler,
                         ψ::mt.QubitCAMPS,
                         n::Int,
                         move_ortho_right::Bool)
  # Search on a copy, so that the real MPS (and its bond dimensions) is untouched when no gate
  # is applied, as in the old code.
  trial = copy(mt.get_mps(ψ))
  orthogonalize!(trial, n)
  id_value = bond_value(de, trial, _ID_2Q, n)
  values = [bond_value(de, trial, M, n) for M in de.matrices]
  min_value, best = findmin(values) # First minimum, as the old code
  min_value < id_value || return 0, 0.0, 0

  mps = mt.get_mps(ψ)
  orthogonalize!(mps, n) # The old applyGateNN! always re-centres at n first
  mt.transform!(ψ, de.gates[best], [n, n + 1];
                move_orthogonality_right = move_ortho_right, de.svd_kwargs...)
  return 1, id_value - min_value, best
end

# == Sweeps =============================================================================

"""```SnakeDisentangler(local_disentangler)```

The old `:snake` sweep: bonds 1, …, N−1, then N−1, …, 1, skipping bonds whose dimension was 1
in a snapshot taken at the start of each direction. Gates are applied with the orthogonality
centre on the left site (old `leftOC = true`). Returns `(number of gates, total gain, info)`,
with `info` the `(gate index, bond)` of each applied gate."""
struct SnakeDisentangler{L <: mt.AbstractLocalDisentangler} <: mt.AbstractStandAloneDisentangler
  local_disentangler::L
end

"""```BrickworkDisentangler(local_disentangler)```

The old `:brickwork` sweep: bonds 1, 3, 5, …, then the even bonds in reverse, with no
bond-dimension skip. Same return values as `SnakeDisentangler`."""
struct BrickworkDisentangler{L <: mt.AbstractLocalDisentangler} <: mt.AbstractStandAloneDisentangler
  local_disentangler::L
end

function sweep_bonds!(de::mt.AbstractLocalDisentangler, ψ::mt.QubitCAMPS, bonds, info)
  num, gain = 0, 0.0
  for b in bonds
    n_b, gain_b, gate = mt.disentangle!(de, ψ, b, false)
    n_b == 0 && continue
    num += n_b
    gain += gain_b
    push!(info, (gate, b))
  end
  return num, gain
end

function mt.disentangle!(de::SnakeDisentangler, ψ::mt.QubitCAMPS)
  N = mt.nsites(ψ)
  info = Tuple{Int, Int}[]
  num, gain = 0, 0.0
  for bonds in (1:N-1, N-1:-1:1)
    dims = linkdims(mt.get_mps(ψ)) # Snapshot per direction, as the old code
    n_d, gain_d = sweep_bonds!(de.local_disentangler, ψ, [b for b in bonds if dims[b] != 1], info)
    num += n_d
    gain += gain_d
  end
  return num, gain, info
end

function mt.disentangle!(de::BrickworkDisentangler, ψ::mt.QubitCAMPS)
  N = mt.nsites(ψ)
  info = Tuple{Int, Int}[]
  num, gain = sweep_bonds!(de.local_disentangler, ψ, [1:2:N-1; reverse(2:2:N-1)], info)
  return num, gain, info
end

"""```MinDiffDisentangler(disentangler, max_iter; min_diff = 0.0)```

Repeat a sweep up to `max_iter` times, stopping after the first sweep whose total gain is
≤ `min_diff` (the old CliffordMPS rule; upstream's `IterativeDisentangler` stops on gain
< `min_gain`). Returns `(number of gates, total gain, per-sweep results)`."""
struct MinDiffDisentangler{D <: mt.AbstractStandAloneDisentangler} <: mt.AbstractStandAloneDisentangler
  disentangler::D
  max_iter::Int
  min_diff::Float64
end

MinDiffDisentangler(de::mt.AbstractStandAloneDisentangler, max_iter::Integer; min_diff::Real = 0.0) =
  MinDiffDisentangler(de, Int(max_iter), Float64(min_diff))

function mt.disentangle!(de::MinDiffDisentangler, ψ::mt.QubitCAMPS)
  sweeps = Tuple{Int, Float64, Any}[]
  for _ in 1:de.max_iter
    num, gain, info = mt.disentangle!(de.disentangler, ψ)
    push!(sweeps, (num, gain, info))
    gain <= de.min_diff && break
  end
  return sum(first, sweeps; init = 0), sum(s -> s[2], sweeps; init = 0.0), sweeps
end

"""```sweep_disentangler(criterion, strategy, max_iter; min_diff = 0.0, svd_kwargs...)```

The MagicTensors equivalent of the old CliffordMPS call
`disentangle!(camps, DisentangleStrategy(strategy), max_iter;
criterion = DisentangleCriterion(criterion), min_diff, svd_kwargs...)`, as a
`MagicTensors.AbstractStandAloneDisentangler`: apply it with `MagicTensors.disentangle!(de, ψ)`.
`criterion` is `:chi3` or `:entangle`, `strategy` is `:snake` or `:brickwork` (the old `:full`,
`:radius` and `:none` are not ported)."""
function sweep_disentangler(criterion::Symbol, strategy::Symbol, max_iter::Integer;
                            min_diff::Real = 0.0, svd_kwargs...)
  bond = BondSearchDisentangler(criterion; svd_kwargs...)
  sweep = if strategy == :snake
    SnakeDisentangler(bond)
  elseif strategy == :brickwork
    BrickworkDisentangler(bond)
  else
    throw(ArgumentError("strategy must be :snake or :brickwork, got :$strategy"))
  end
  return MinDiffDisentangler(sweep, max_iter; min_diff = min_diff)
end
