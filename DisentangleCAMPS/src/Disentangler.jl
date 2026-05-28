using QuantumClifford

"After mapping the Pauli string P through the Clifford operator C (ie. P' = C^†PC),
determine the type of action the result has on the state |m⟩^(⊗k) |0⟩^(⊗(N-k)):
- ```:disentanglable```, trivial action on the first k qubits and\
nontrivial on at least one other: action of a phase gate e^iϕP can be\
incorporated without increasing the MPS bond dimension
- ```:trivial```, trivial action on all qubits: action of phase gate can be discarded
- ```:logical```, nontrivial action on one of the first k qubits: MPS bond dimension must increase."
function paulinature(k::Integer, C::CliffordOperator, P::PauliOperator)
  n = length(P)
  Q = apply(P, C)
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

"Map a Pauli string P through a Clifford operator C, ie. compute C^†PC.\
Warning: not an in-place method!"
function apply(P::PauliOperator, C::CliffordOperator)
  P_tableau = QuantumClifford.Tableau([P])
  apply!(Stabilizer(P_tableau), inv(C))
  return P_tableau[1]
end

"Given a Pauli string P, a Clifford operator C and the number k of free qubits,\
build an analytical disentangling Clifford circuit."
function disentangler(k::Integer, C::CliffordOperator, P::PauliOperator)
  Q = apply(P, C)
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

disentangler(k::Integer, P::PauliOperator) = disentangler(k, one(CliffordOperator, length(P)), P)

"Routine for ```disentangler``` to use if the first free qubit with nontrivial action is not number k+1"
function swapqubits(P::PauliOperator, i::Integer, j::Integer)
  swapij = one(CliffordOperator, length(P))
  apply!(swapij, tSWAP, [i,j])
  return apply(P, swapij), swapij
end

"Routine for ```disentangler``` to use if the first free qubit with nontrivial action is acted on by a Y"
function reducetoX(P::PauliOperator, i::Integer)
  phasei = one(CliffordOperator, length(P))
  apply!(phasei, tPhase, [i]) # Use 3 phase gates s.t. Y -> X instead of -X
  return apply(P, inv(phasei)), phasei
end

findfirstfreeXY(P::PauliOperator, k::Integer) = findfirst(xbit(P)[k+1:end]) + k

"Routine for ```disentangler``` to use to build the disentangling circuit from Liu and Clark (2025)."
function build_D(Q::PauliOperator, i::Integer)
  tCX = tCNOT
  tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
  tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

  D = one(CliffordOperator, length(Q))

  for j in eachindex(Q)
    if j != i
      if Q[j] == (true, false) # X
        apply!(D, tCX, [i, j])
      elseif Q[j] == (true, true) # Y
        apply!(D, tCY, [i, j])
      elseif Q[j] == (false, true) # Z
        apply!(D, tCZ, [i, j])
      end
    end
  end
  # println("")
  return D
end