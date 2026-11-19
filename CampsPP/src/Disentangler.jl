using QuantumClifford

"""```paulinature(k, C, P)```

Classify the mapped string Q = C†PC (C the CAMPS frame, k the number of magic qubits, \
which are qubits 1:k; qubits k+1:N are free, in |0⟩ in the MPS) by how exp(iϕQ) acts on \
the MPS:
- ```:disentanglable```: Q has X or Y on at least one free qubit. exp(iϕP) can be absorbed \
  as a Clifford frame update plus a new magic qubit k+1, with no growth of the bond \
  dimension (see ```disentangler```).
- ```:trivial```: identity on the magic qubits and only I or Z on the free ones, which are \
  in |0⟩, so exp(iϕQ) is a global phase and can be dropped.
- ```:logical```: X, Y or Z on some magic qubit and only I or Z on the free ones: \
  exp(iϕQ) must be applied to the MPS, which can increase the bond dimension."""
paulinature(k::Integer, C::CliffordOperator, P::PauliOperator) = paulinature(k, apply(P, C))

"```paulinature(k, Q)```

As `paulinature(k, C, P)`, for an already mapped string Q = C^†PC (identity frame)."
function paulinature(k::Integer, Q::PauliOperator)
  n = length(Q)
  XQ = xbit(Q)
  ZQ = zbit(Q)
  if XQ[k+1:n] != zeros(n-k)
    return :disentanglable
  elseif ZQ[1:k] == zeros(k) && XQ[1:k] == zeros(k)
    return :trivial
  else
    return :logical
  end
end

"""```apply(P, C)```

The Pauli string P mapped through the Clifford operator C, i.e. C†PC, with its sign. \
Returns a new string; P and C are not modified."""
function apply(P::PauliOperator, C::CliffordOperator)
  P_tableau = QuantumClifford.Tableau([P])
  apply!(Stabilizer(P_tableau), inv(C))
  return P_tableau[1]
end

"""```_apply_on_sites!(C, op, sites)```

Compose the Clifford `op`, acting on the qubits `sites`, after C (C ← op·C, phases tracked), in
place. QuantumClifford has no `apply!(::CliffordOperator, ::CliffordOperator, sites)`; this
replaces the one the old vendored CliffordMPS defined (plan, design correction 3).
`Stabilizer(C.tab)` wraps C's tableau without copying it, so the write goes to C."""
function _apply_on_sites!(C::CliffordOperator, op::CliffordOperator, sites::AbstractVector{<:Integer})
  apply!(Stabilizer(C.tab), sites, op)
  return C
end

"""```disentangler(k, C, P)```

The analytical disentangling Clifford of Liu and Clark (2025) for a rotation exp(iϕP) on \
a CAMPS with frame C and k magic qubits, where Q = C†PC is ```:disentanglable``` (see \
```paulinature```). Returns ```(D, s)``` with s = ±1 the sign of Q: the rotation is then \
applied as the frame update C ← C·D and exp(i·sϕ·X_{k+1}) on qubit k+1 of the MPS (as in \
```apply!(ψ, k, P, ϕ)```), since D†QD = s·X_{k+1} (test/CONVENTIONS.md §4). D is built from at most one phase gate (if Q has Y on the \
chosen qubit), one SWAP (moving the first free qubit with X or Y to k+1) and controlled \
Paulis from qubit k+1 (```build_D```)."""
disentangler(k::Integer, C::CliffordOperator, P::PauliOperator) = disentangler(k, apply(P, C))

"""```disentangler(k, Q)```

As ```disentangler(k, C, P)```, for an already mapped string Q = C†PC (identity frame)."""
function disentangler(k::Integer, Q::PauliOperator)
  sign = (Q.phase[1] == 0x0) ? +1 : -1
  i = findfirstfreeXY(Q, k)
  is_y = zbit(Q)[i]
  Dtot = one(CliffordOperator, length(Q))
  
  if is_y
    Q, phase = reducetoX(Q, i)
    apply_right!(Dtot, phase)
  end
  
  if i != k+1
    Q, swap = swapqubits(Q, i, k+1)
    apply_right!(Dtot, swap)
  end
  
  Dmain = build_D(Q, k+1)
  apply_right!(Dtot, Dmain)

  return Dtot, sign
end

"""```swapqubits(P, i, j)```

For ```disentangler```, when the first free qubit with X or Y is not k+1: the SWAP of qubits \
i and j as a Clifford, and P mapped through it. Returns ```(mapped P, SWAP)```."""
function swapqubits(P::PauliOperator, i::Integer, j::Integer)
  swapij = one(CliffordOperator, length(P))
  _apply_on_sites!(swapij, tSWAP, [i,j])
  return apply(P, swapij), swapij
end

"""```reducetoX(P, i)```

For ```disentangler```, when P has Y on qubit i: the phase gate S = diag(1, i) on qubit i, \
and S·P·S† (```apply(P, S†)```), which turns the Y into X with the overall sign flipped \
(S·Y·S† = −X; checked against dense matrices). Returns ```(S·P·S†, S)```."""
function reducetoX(P::PauliOperator, i::Integer)
  phasei = one(CliffordOperator, length(P))
  _apply_on_sites!(phasei, tPhase, [i]) # Returned string is S·P·S†: +Y → −X
  return apply(P, inv(phasei)), phasei
end

"The first free qubit (index > k) on which P has X or Y; errors if there is none."
findfirstfreeXY(P::PauliOperator, k::Integer) = findfirst(xbit(P)[k+1:end]) + k

"""```build_D(Q, i)```

For ```disentangler```: the disentangling circuit of Liu and Clark (2025) for a string Q \
with X on qubit i, a controlled-X, -Y or -Z from qubit i onto every other qubit where Q \
has X, Y or Z."""
function build_D(Q::PauliOperator, i::Integer)
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

  D = one(CliffordOperator, length(Q))

  for j in eachindex(Q)
    if j != i
      if Q[j] == (true, false) # X
        _apply_on_sites!(D, tCX, [i, j])
      elseif Q[j] == (true, true) # Y
        _apply_on_sites!(D, tCY, [i, j])
      elseif Q[j] == (false, true) # Z
        _apply_on_sites!(D, tCZ, [i, j])
      end
    end
  end
  return D
end