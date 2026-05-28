using ITensors
using Random

@doc """
    sEntropy(psi::MPS; α=2.0)

Brute force computation of stabilizer α-Renyi Entropy [1].
For α=-1.0 cumputes the stabilizer linear entropy.

[1] ... Leone, Oliviero, Hamma; PRL (2022); DOI: 10.1103/PhysRevLett.128.050402
"""
function sEntropy(psi::MPS; α=2.0, filtered=false)
    return sEntropy(psi, [α]; filtered=filtered)[1]
end

function sEntropy(psi::MPS, α::Vector; filtered=false)
    paulibasis = ["I" "X" "Y" "Z"]
    sites = siteinds(psi)
    N = length(sites)
    d = 2^N
    sumVal = 1.0

    paulistrings = [""]
    for n in 1:N
        paulistrings = kron(paulistrings, paulibasis)
    end

    if filtered
        paulistrings = paulistrings[2:end]
        sumVal = 1.0 - 1/d
    end


    Ξ = []
    for paulistring in paulistrings
        paulioperator = op(paulistring[1:1],sites[1])
        for n in 2:N
            paulioperator *= op(paulistring[n:n],sites[n])
        end
        val = real(inner(psi',MPO(paulioperator, siteinds(psi)),psi))
        push!(Ξ, val^2 / d)
    end

    if !(sum(Ξ) ≈ sumVal)
        @warn "Sum of paulisting expectationation should be close to $(sumVal), but it is $(sum(Ξ))!"
    end

    results = []

    Ξ = Ξ[findall(>(0.0),Ξ)]
    for a in α
        if a == 1.0
            push!(results,-sum(Ξ .* log2.(Ξ)) - N)
        elseif a == -1.0
            push!(results,1 - 2^N * sum(Ξ.^2))
        else
            push!(results,1/(1-a)*log2(sum(Ξ.^a)) - N)
        end
    end
    return results
end

@doc """
    sEntropy(psi::MPS, samples::Int; α=[2.0])

Compute the stabilizer α-Renyi entropy by sampling the Paulistring distribution [1].
For α=-1.0 cumputes the stabilizer linear entropy.

[1] ... Lami, Collura; arXiv:2303.05536; DOI: 10.48550/arXiv.2303.05536
"""
function sEntropy(psi::MPS, samples::Integer, α::Vector; filtered=false)
    psi = orthogonalize(psi,1)
    normalize!(psi)
    s = siteinds(psi)
    N = length(s)

    Π = []
    for i in 1:samples
        push!(Π, get_pauli_sample(psi::MPS; filtered=filtered))
    end

    results = []

    Π = Π[findall(>(0.0),Π)]
    for a in α
        if a == 1.0
            push!(results,-sum(log2.(Π))/samples - N)
        elseif a == -1.0
            push!(results,1-2^N*(sum(Π)/samples))
        else
            push!(results, 1/(1-a)*log2(sum(Π.^(a-1))/samples) - N)
        end
    end

    return results
end

function sEntropy(psi::MPS, samples::Integer; α=2.0, filtered=false)
    return sEntropy(psi, samples, [α]; filtered=filtered)[1]
end

@doc """
    get_pauli_sample(psi::MPS)

Sample a pauli string [1]. The MPS psi has to be normalized and the
orthogonality has to be at site 1.

    orthogonalize!(psi,1)
    normalize!(psi)

[1] ... Lami, Collura; arXiv:2303.05536; DOI: 10.48550/arXiv.2303.05536
"""
function get_pauli_sample(psi::MPS; filtered=false)
    s = siteinds(psi)
    N = length(s)

    paulis = ["I", "X", "Y", "Z"]
    L = ITensor(ComplexF64, 1.0)
    Π = 1.0
    fullIdentity = true

    for n ∈ 1:N
        r = rand()
        for pauli in paulis
            l = L
            l = l * psi[n] * op(pauli,s[n]) * prime(dag(psi[n]))
            p = real(Array(l * dag(l))[1,1]) / 2.0
            r -= p
            if r ≤ 0.0
                if fullIdentity && pauli != "I"
                    fullIdentity = false
                end
                L = l / √(2.0*p)
                Π *= p
                break
            end
        end
    end
    if filtered && fullIdentity
        return get_pauli_sample(psi; filtered=true)
    end
    return Π
end
