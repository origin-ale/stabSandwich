# The kicked Ising circuit of test/CONVENTIONS.md §2 "Kicked Ising", built from the spec
# alone (not from `kicked_ising_circuit`), and its dense ⟨Z_j⟩ profile. Shared by the CampsPP
# tests (through helpers.jl) and test/TestKimScript.jl. Needs `DR` (the dense reference
# module, CampsPP/test/DisentangleCAMPS/DenseReference.jl) in the including module.
#
# Spec, in the exp(iϕP) convention: per layer, (ZZ, −π/4) on every bond (the bricklayer
# bonds (1,2), (3,4), … then (2,3), (4,5), …, then the wrap bond (1, N) for N ≥ 3 only), then
# (Z_j, −(π/2 + ε_j)/2) on j = 1…N, then (Y_j, −π/4) on j = 1…N. Initial state Néel
# |0101…⟩: qubit j is |1⟩ for even j.

"Bonds of one KIM layer, in gate order; `wrap = false` drops the wrap bond (an open chain)."
function kim_spec_bonds(N; wrap = true)
  bonds = vcat([[i, i + 1] for i in 1:2:N-1], [[i, i + 1] for i in 2:2:N-1])
  return wrap && N >= 3 ? [bonds; [[1, N]]] : bonds
end

"Gates per layer, from the spec formula: L = B + 2N, B = N for N ≥ 3, N − 1 for N ≤ 2."
kim_spec_L(N) = (N >= 3 ? N : N - 1) + 2N

"The Néel state's `onebitinds`: the even qubits."
kim_neel(N) = collect(2:2:N)

"N-character string with `sym` on the qubits `qs` and I elsewhere."
kim_pstring(N, sym::Char, qs) = String([q in qs ? sym : 'I' for q in 1:N])

"""
`kim_spec_circuit(ε, t, N; wrap = true)`: t KIM layers as (strings, phases), strings over
IXYZ (character j = qubit j), phases in the exp(iϕP) convention. ε is a number (every qubit)
or a length-N vector (ε_j on qubit j).
"""
function kim_spec_circuit(εs::AbstractVector, t, N; wrap = true)
  @assert length(εs) == N
  bonds = kim_spec_bonds(N; wrap)
  strings = [[kim_pstring(N, 'Z', b) for b in bonds];
             [kim_pstring(N, 'Z', (j,)) for j in 1:N];
             [kim_pstring(N, 'Y', (j,)) for j in 1:N]]
  phases = [fill(-π/4, length(bonds)); [-(π/2 + ε) / 2 for ε in εs]; fill(-π/4, N)]
  return repeat(strings, t), repeat(Float64.(phases), t)
end
kim_spec_circuit(ε::Real, t, N; kw...) = kim_spec_circuit(fill(ε, N), t, N; kw...)

"""
`kim_dense_profile(strings, phases, N, layer; ones = kim_neel(N))`: dense ⟨Z_j⟩ for every
qubit j, at t = 0 and after every `layer` gates: a (t + 1) × N matrix, row l + 1 = after
layer l. `layer` is the number of gates per layer (it differs from kim_spec_L(N) for an
open chain).
"""
function kim_dense_profile(strings, phases, N, layer; ones = kim_neel(N))
  ψ = DR.basis_state(N, ones)
  zs() = [DR.expect(ψ, kim_pstring(N, 'Z', (j,))) for j in 1:N]
  rows = [zs()]
  for (i, (P, ϕ)) in enumerate(zip(strings, phases))
    DR.rotate!(ψ, P, ϕ)
    i % layer == 0 && push!(rows, zs())
  end
  return permutedims(reduce(hcat, rows))
end

"Dense ⟨Z_j⟩ profile of t spec layers with deviation(s) ε from the Néel state ((t + 1) × N)."
kim_dense_profile(ε, t, N; wrap = true, ones = kim_neel(N)) =
  kim_dense_profile(kim_spec_circuit(ε, t, N; wrap)..., N,
                    length(first(kim_spec_circuit(ε, 1, N; wrap))); ones)
