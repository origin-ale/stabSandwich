# Stabilizer-tableau reference for circuits whose angles are all multiples of π/4, at any N
# (T6, large-N Clifford tests). Written from test/CONVENTIONS.md §2 only, not from the code:
#
#   for U = exp(iϕP) and a Pauli string A,  U A U† = A                  if A commutes with P,
#                                            U A U† = exp(2iϕP)·A
#                                                   = cos 2ϕ·A + sin 2ϕ·iPA  if A anticommutes.
#
# At a multiple of π/4 one of cos 2ϕ, sin 2ϕ is 0 and the other ±1, so each stabilizer
# generator maps to one signed Pauli string. The state is stored as its N generators.
#
# QuantumClifford is used only for Pauli products (`*`, phases tracked) and for `expect`
# on a stabilizer state. Both are checked against the dense reference by the self-test in
# CampsPP/test/DisentangleCAMPS/TestCliffordLargeN.jl (all 4^N strings at N ≤ 3).

module StabilizerReference

using QuantumClifford

export anticommutes, basis_generators, rotate_generators!, run_generators, stab_expect

"Do the Pauli strings P and A anticommute? (Symplectic product Σⱼ xᴾⱼ zᴬⱼ + zᴾⱼ xᴬⱼ odd.)"
anticommutes(P::PauliOperator, A::PauliOperator) =
  isodd(count((xbit(P) .& zbit(A)) .⊻ (zbit(P) .& xbit(A))))

"Generators of the basis state with `ones` in |1⟩ (CONVENTIONS.md §1): −Zⱼ for j in ones, +Zⱼ otherwise."
function basis_generators(N::Integer, ones)
  map(1:N) do j
    Z = PauliOperator(0x0, falses(N), [i == j for i in 1:N])
    j in ones ? -Z : Z
  end
end

"""
Replace each generator A by U A U† for U = exp(iϕP), ϕ a multiple of π/4 (see the header).
P may carry a sign: exp(iϕ(−P)) = exp(−iϕP) follows from the same formula.
"""
function rotate_generators!(gens::Vector{<:PauliOperator}, P::PauliOperator, ϕ::Real)
  c, s = cos(2ϕ), sin(2ϕ)
  ci, si = round(Int, c), round(Int, s)
  (abs(c - ci) < 1e-9 && abs(s - si) < 1e-9) ||
    throw(ArgumentError("ϕ = $ϕ is not a multiple of π/4"))
  for (j, A) in enumerate(gens)
    anticommutes(P, A) || continue
    B = si == 0 ? A : im * (P * A)              # cos 2ϕ·A or i·P·A
    B = (ci + si) == 1 ? B : -B                 # the nonzero one of cos 2ϕ, sin 2ϕ is ±1
    B.phase[] in (0x0, 0x2) || error("non-Hermitian generator: bug in the reference")
    gens[j] = B
  end
  return gens
end

"Copy of `gens` after the circuit (paulis, phis), gates[1] first (CONVENTIONS.md §2)."
function run_generators(gens, paulis, phis)
  gens = copy(gens)
  for (P, ϕ) in zip(paulis, phis)
    rotate_generators!(gens, P, ϕ)
  end
  return gens
end

"⟨P⟩ on the stabilizer state with generators `gens` (−1, 0 or 1)."
stab_expect(gens, P::PauliOperator) = stab_expect(Stabilizer(gens), P)
function stab_expect(S::Stabilizer, P::PauliOperator)
  ev = QuantumClifford.expect(P, S)
  isreal(ev) || error("non-Hermitian observable")
  return real(ev)
end

end # module
