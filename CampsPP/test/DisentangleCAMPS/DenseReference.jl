"""
    DenseReference

Exact state-vector simulator for a few qubits (N ≲ 12), written independently of the
package it is used to check. Conventions follow test/CONVENTIONS.md §1–2:
qubits 1…N, |0⟩ the +1 eigenstate of Z, Y = [[0, -i], [i, 0]], a gate (P, ϕ) is
exp(iϕP) = cos ϕ + i sin ϕ · P, and the first gate of a list acts first.

Indexing: the amplitude of |b₁ b₂ … b_N⟩ (bⱼ ∈ {0,1} the state of qubit j) is stored at
`ψ[1 + Σⱼ bⱼ 2^(j-1)]`, i.e. qubit 1 is the least significant bit. This is the
column-major layout of an array ψ[b₁,…,b_N] with qubit 1 as the first index, and the
matrix of A₁ ⊗ … ⊗ A_N (Aⱼ on qubit j) is `kron(A_N, …, A_1)`.
"""
module DenseReference

export zero_state, basis_state, apply_pauli, rotate, rotate!, run_circuit,
    expect, all_paulis, all_expectations, bloch

"Number of qubits of a state vector (its length must be a power of two)."
function nqubits(ψ::AbstractVector)
    N = trailing_zeros(length(ψ))
    length(ψ) == 1 << N || throw(ArgumentError("state length $(length(ψ)) is not a power of 2"))
    return N
end

"""
    parse_pauli(pauli, N) -> (sign, xmask, zmask, ny)

Bit j-1 of `xmask` is set if character j is X or Y (the bit is flipped); bit j-1 of
`zmask` is set if it is Z or Y (a factor (-1)^bⱼ); `ny` counts the Y's, each of which
adds a factor i since Y = iXZ. `sign` is ±1 from an optional leading '+' or '-'.
"""
function parse_pauli(pauli::AbstractString, N::Integer)
    chars = collect(pauli)
    sign = 1
    if !isempty(chars) && chars[1] in ('+', '-')
        sign = chars[1] == '-' ? -1 : 1
        popfirst!(chars)
    end
    length(chars) == N ||
        throw(ArgumentError("Pauli string \"$pauli\" has $(length(chars)) letters, expected $N"))
    xmask = 0; zmask = 0; ny = 0
    for (j, c) in enumerate(chars)
        bit = 1 << (j - 1)
        if c == 'X'
            xmask |= bit
        elseif c == 'Z'
            zmask |= bit
        elseif c == 'Y'
            xmask |= bit; zmask |= bit; ny += 1
        elseif c != 'I'
            throw(ArgumentError("invalid character '$c' in Pauli string \"$pauli\""))
        end
    end
    return sign, xmask, zmask, ny
end

"`zero_state(N)`: |0…0⟩ on N qubits."
zero_state(N::Integer) = basis_state(N, Int[])

"`basis_state(N, ones)`: the computational basis state with the qubits in `ones` in |1⟩."
function basis_state(N::Integer, ones)
    idx = 0
    for j in ones
        1 <= j <= N || throw(ArgumentError("qubit $j is outside 1:$N"))
        idx |= 1 << (j - 1)
    end
    ψ = zeros(ComplexF64, 1 << N)
    ψ[idx + 1] = 1
    return ψ
end

"""
    apply_pauli(ψ, pauli) -> P|ψ⟩

On a basis state: P|b⟩ = sign · i^ny · (-1)^popcount(b & zmask) · |b ⊻ xmask⟩.
Per qubit this is X|b⟩ = |1-b⟩, Z|b⟩ = (-1)^b|b⟩, Y|0⟩ = i|1⟩, Y|1⟩ = -i|0⟩.
"""
function apply_pauli(ψ::AbstractVector, pauli::AbstractString)
    sign, xmask, zmask, ny = parse_pauli(pauli, nqubits(ψ))
    pref = sign * (im)^ny
    out = zeros(ComplexF64, length(ψ))
    for b in 0:length(ψ)-1
        phase = isodd(count_ones(b & zmask)) ? -pref : pref
        out[(b ⊻ xmask) + 1] = phase * ψ[b + 1]
    end
    return out
end

"`rotate!(ψ, pauli, ϕ)`: ψ ← exp(iϕP)ψ = cos ϕ · ψ + i sin ϕ · Pψ, in place."
function rotate!(ψ::AbstractVector, pauli::AbstractString, ϕ::Real)
    Pψ = apply_pauli(ψ, pauli)
    c, s = cos(ϕ), sin(ϕ)
    @. ψ = c * ψ + im * s * Pψ
    return ψ
end

"`rotate(ψ, pauli, ϕ)`: exp(iϕP)|ψ⟩ as a new vector."
rotate(ψ::AbstractVector, pauli::AbstractString, ϕ::Real) =
    rotate!(Vector{ComplexF64}(ψ), pauli, ϕ)

"`run_circuit(ψ, paulis, phis)`: apply the gates (paulis[t], phis[t]) for t = 1, 2, … in that order."
function run_circuit(ψ::AbstractVector, paulis, phis)
    length(paulis) == length(phis) ||
        throw(ArgumentError("$(length(paulis)) Pauli strings but $(length(phis)) angles"))
    out = Vector{ComplexF64}(ψ)
    for (P, ϕ) in zip(paulis, phis)
        rotate!(out, P, ϕ)
    end
    return out
end

"`expect(ψ, pauli)`: ⟨ψ|P|ψ⟩ as a real number, without renormalising ψ."
function expect(ψ::AbstractVector, pauli::AbstractString)
    Pψ = apply_pauli(ψ, pauli)
    val = sum(conj(ψ[k]) * Pψ[k] for k in eachindex(ψ))
    abs(imag(val)) <= 1e-9 || error("⟨$pauli⟩ = $val has a non-negligible imaginary part")
    return real(val)
end

"""
    all_paulis(N)

The 4^N unsigned strings. Entry k (1-based) has, on qubit j, the letter "IXYZ"[d+1]
with d the j-th base-4 digit of k-1 (qubit 1 = least significant digit, so it varies
fastest): for N = 2 the order is II, XI, YI, ZI, IX, XX, …, ZZ.
"""
all_paulis(N::Integer) =
    [join("IXYZ"[((k >> (2 * (j - 1))) & 3) + 1] for j in 1:N) for k in 0:4^N-1]

"`all_expectations(ψ)`: ⟨ψ|P|ψ⟩ for every P in `all_paulis(N)`, in the same order."
all_expectations(ψ::AbstractVector) = [expect(ψ, P) for P in all_paulis(nqubits(ψ))]

"`bloch(ψ, j)`: (⟨X_j⟩, ⟨Y_j⟩, ⟨Z_j⟩)."
function bloch(ψ::AbstractVector, j::Integer)
    N = nqubits(ψ)
    1 <= j <= N || throw(ArgumentError("qubit $j is outside 1:$N"))
    single(c) = String([i == j ? c : 'I' for i in 1:N])
    return (expect(ψ, single('X')), expect(ψ, single('Y')), expect(ψ, single('Z')))
end

end # module
