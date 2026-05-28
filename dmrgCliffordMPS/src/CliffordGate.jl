using CliffordMPS
using QuantumClifford
using ITensors
using StatsBase

import Base

# -- Abstract Clifford Indices Sets ---------------------------------------------------------------

abstract type CliffordIndexSets end

function _check_valid_num(cliff_set::CliffordIndexSets, num)
    @assert(num>=0 && num<length(cliff_set))
end

function _mod_div(num::Int, cliff_set::CliffordIndexSets)
    n = num % length(cliff_set)
    new_num = div(num,length(cliff_set))
    ind = clifford_type(cliff_set)(cliff_set, n)
    return new_num, ind
end

# -- 1Qubit Clifford Gate Numbering ---------------------------------------------------------------

abstract type Clifford1IndexSets <: CliffordIndexSets end
clifford_type(::Clifford1IndexSets) = Clifford1Index

struct Clifford1IndexSet{x} <: Clifford1IndexSets
end

Clifford1IndexSet(x) = Clifford1IndexSet{x}()

"""
Structured index for a one-qubit Clifford gate.

The fields `a` and `b` encode the decomposition used by this package's
enumeration scheme. Use `clifford_gate(ind)` to obtain the corresponding gate.
"""
struct Clifford1Index
    a::Int
    b::Int
    function Clifford1Index(a::Int, b::Int)
        _check_valid_index(a, b)
        return new(a, b)
    end
end

function _check_valid_index(a::Int, b::Int)
    @assert(a>=0 && a<6)
    @assert(b>=0 && b<4)
    return true
end

function Clifford1Index(num::Int)
    return Clifford1Index(Clifford1IndexSet(:all), num)
end
function clifford_number(ind::Clifford1Index)
    return clifford_number(Clifford1IndexSet(:all), ind)
end


Base.length(::Clifford1IndexSet{:all}) = 6*4
function Clifford1Index(::Clifford1IndexSet{:all}, num::Int)
    b = num % 4
    num = div(num, 4)
    a = num % 6
    num = div(num, 6)
    @assert(num == 0)
    return Clifford1Index(a, b)
end
function clifford_number(::Clifford1IndexSet{:all}, ind::Clifford1Index)
    return ind.a*4+ind.b
end

_rotations = [Clifford1Index(0,0), Clifford1Index(3,3), Clifford1Index(4,0)]
Base.length(::Clifford1IndexSet{:rotation}) = 3
function Clifford1Index(cliff_set::Clifford1IndexSet{:rotation}, num::Int)
    _check_valid_num(cliff_set, num)
    return _rotations[num+1]
end
function clifford_number(::Clifford1IndexSet{:rotation}, ind::Clifford1Index)
    num = nothing
    for i in 1:length(_rotations)
        if _rotations[i] == ind
            num = i-1
            break
        end
    end
    if num === nothing
        throw(ErrorException("Clifford1Index is not a rotation."))
    end
    return num
end

Base.length(::Clifford1IndexSet{:a}) = 6
function Clifford1Index(::Clifford1IndexSet{:a}, num::Int)
    a = num % 6
    num = div(num, 6)
    @assert(num == 0)
    return Clifford1Index(a, 0)
end
function clifford_number(::Clifford1IndexSet{:a}, ind::Clifford1Index)
    @assert(ind.b == 0)
    return ind.a
end

# -- 2Qubit Clifford Gate Numbering ---------------------------------------------------------------

abstract type Clifford2IndexSets <: CliffordIndexSets end
clifford_type(::Clifford2IndexSets) = Clifford2Index

struct Clifford2IndexSet{x} <: Clifford2IndexSets
end

Clifford2IndexSet(x) = Clifford2IndexSet{x}()


"""
Structured index for a two-qubit Clifford gate.

The representation separates the local Clifford factors from the entangling
core so the package can enumerate different gate families consistently.
"""
struct Clifford2Index
    r::Int
    c1::Clifford1Index
    c2::Clifford1Index
    s1::Clifford1Index
    s2::Clifford1Index
end

function Clifford2Index(r::Int, a1::Int, b1::Int, a2::Int, b2::Int, s1::Int, s2::Int)
    _check_valid_index(r, a1, b1, a2, b2, s1, s2)
    return Clifford2Index(r, Clifford1Index(a1, b1), Clifford1Index(a2, b2),
        Clifford1Index(Clifford1IndexSet(:rotation),s1),
        Clifford1Index(Clifford1IndexSet(:rotation),s2))
end

function _check_valid_index(r::Int, a1::Int, b1::Int, a2::Int, b2::Int, s1::Int, s2::Int)
    @assert(r>=0 && r<4)
    @assert(a1>=0 && a1<6)
    @assert(b1>=0 && b1<4)
    @assert(a2>=0 && a2<6)
    @assert(b2>=0 && b2<4)
    @assert(s1>=0 && s1<3)
    @assert(s2>=0 && s2<3)
    if r<2 && (s1!=0 || s2!=0)
        throw(ArgumentError("s1 and s2 must be 0 for r<2!"))
    end
    return true
end

function _index_numbers(ind::Clifford2Index)
    r = ind.r
    a1 = ind.c1.a
    b1 = ind.c1.b
    a2 = ind.c2.a
    b2 = ind.c2.b
    s1 = clifford_number(Clifford1IndexSet(:rotation), ind.s1)
    s2 = clifford_number(Clifford1IndexSet(:rotation), ind.s2)
    return (r, a1, b1, a2, b2, s1, s2)
end


function Clifford2Index(num::Int)
    return Clifford2Index(Clifford2IndexSet(:all), num)
end

function clifford_number(ind::Clifford2Index)
    return clifford_number(Clifford2IndexSet(:all), ind)
end


Base.length(::Clifford2IndexSet{:all}) = 2*6*4*6*4 + 2*6*4*6*4*3*3
function Clifford2Index(::Clifford2IndexSet{:all}, num::Int)
    short = false
    r = nothing
    c1 = nothing
    c2 = nothing
    s1 = nothing
    s2 = nothing
    if div(num, 6*4*6*4) < 2
        r = div(num, 6*4*6*4)
        num -= r*6*4*6*4
        short = true
    else
        num -= 2*6*4*6*4
        t = div(num, 6*4*6*4*3*3)
        num -= t*6*4*6*4*3*3
        r = t + 2
        short = false
    end
    if short
        s2 = Clifford1Index(Clifford1IndexSet(:rotation),0)
        s1 = Clifford1Index(Clifford1IndexSet(:rotation),0)
    else
        num, s2 = _mod_div(num, Clifford1IndexSet(:rotation))
        num, s1 = _mod_div(num, Clifford1IndexSet(:rotation))
    end
    num, c2 = _mod_div(num, Clifford1IndexSet(:all))
    num, c1 = _mod_div(num, Clifford1IndexSet(:all))
    @assert(num == 0)
    return Clifford2Index(r, c1, c2, s1, s2)
end

function clifford_number(::Clifford2IndexSet{:all}, ind::Clifford2Index)
    N_all = length(Clifford1IndexSet(:all))
    N_rot = length(Clifford1IndexSet(:rotation))
    nc1 = clifford_number(Clifford1IndexSet(:all), ind.c1)
    nc2 = clifford_number(Clifford1IndexSet(:all), ind.c2)
    ns1 = clifford_number(Clifford1IndexSet(:rotation), ind.s1)
    ns2 = clifford_number(Clifford1IndexSet(:rotation), ind.s2)
    if ind.r<2
        return ind.r*N_all^2 + nc1*N_all + nc2
    else
        return 2*N_all^2 + (ind.r-2)*N_all^2*N_rot^2 +
            nc1*N_all*N_rot^2 + nc2*N_rot^2 + ns1*N_rot + ns2
    end
end



Base.length(::Clifford2IndexSet{:A}) = 1 + 2*6*6
function Clifford2Index(::Clifford2IndexSet{:A}, num::Int)
    # num ∈ 0 ... 72
    # r ∈ {[1], 2, 3}
    # a1 ∈ {0, 1, 2, 3, 4, 5}
    # a2 ∈ {0, 1, 2, 3, 4, 5}
    if num == 0
        return Clifford2Index(1,0,0,0,0,0,0)
    else
        num -= 1
    end
    s2 = Clifford1Index(Clifford1IndexSet(:rotation),0)
    s1 = Clifford1Index(Clifford1IndexSet(:rotation),0)
    num, c2 = _mod_div(num, Clifford1IndexSet(:a))
    num, c1 = _mod_div(num, Clifford1IndexSet(:a))
    r = num + 2
    return Clifford2Index(r, c1, c2, s1, s2)
end

function clifford_number(::Clifford2IndexSet{:A}, ind::Clifford2Index)
    if ind == Clifford2Index(1,0,0,0,0,0,0)
        return 0
    end
    num = 1
    N_a = length(Clifford1IndexSet(:a))
    num += (ind.r-2) * N_a^2
    num += clifford_number(Clifford1IndexSet(:a), ind.c1) * N_a
    num += clifford_number(Clifford1IndexSet(:a), ind.c2)
end


entangle_set = [576, 1152, 1188, 1260, 2016, 2052, 2124, 3744, 3780, 3852, 6336, 6372, 6408, 7200, 7236, 7272, 8928, 8964, 9000]

Base.length(::Clifford2IndexSet{:entangle}) = length(entangle_set)
function Clifford2Index(::Clifford2IndexSet{:entangle}, num::Int)
    # num ∈ 0 ... length-1
    return Clifford2Index(entangle_set[num+1])
end

function clifford_number(::Clifford2IndexSet{:entangle}, ind::Clifford2Index)
    full_num = clifford_number(ind)
    return findall(x->(x==full_num),entangle_set)[1]-1
end

# -- Draw Clifford Circuits -----------------------------------------------------------------------
    
strA = ["-----", "--H--", "--S--", "-SH--", "-HS--", "-HSH-", "     "]
strB = ["---", "-X-", "-Y-", "-Z-", "   "]
strR=  [("------","      ","------"),
        ("-SWAP-","  |   ","--o---"),
        ("-iSWP-","  |   ","--o---"),
        ("-CNOT-","  |   ","--o---"),]
strS = ["----", "-R--", "-RR-", "    "]

draw(io::IO, ind::Clifford1Index) = print(io, strB[ind.b+1], strA[ind.a+1])
draw(io::IO, ind::Clifford1Index, ::Nothing) = print(io, strB[end], strA[end])
draw(io::IO, ind::Clifford1Index, ::Clifford1IndexSet{:rotation}, ::Nothing) =  print(io, strS[end])
function draw(io::IO, ind::Clifford1Index, ::Clifford1IndexSet{:rotation})
    s = clifford_number(Clifford1IndexSet(:rotation),ind)
    print(io, strS[s+1])
end

function draw(io::IO, ind::Clifford2Index)
    draw(io,ind.c1)
    print(io,strR[ind.r+1][1])
    draw(io,ind.s1,Clifford1IndexSet(:rotation))
    println("")
    draw(io,Clifford1Index(0),nothing)
    print(io,strR[ind.r+1][2])
    draw(io,Clifford1Index(0),Clifford1IndexSet(:rotation),nothing)
    println("")
    draw(io,ind.c2)
    print(io,strR[ind.r+1][3])
    draw(io,ind.s2,Clifford1IndexSet(:rotation))
    println("")
end

function Base.show(io::IO, ind::Clifford1Index)
    println(io, "Clifford1Index",(ind.a, ind.b))
    draw(io, ind)
end

function Base.show(io::IO, ind::Clifford2Index)
    println(io, "Clifford2Index",_index_numbers(ind))
    draw(io, ind)
end

# -- Stabilizer Clifford Gates --------------------------------------------------------------------

tI = tId1
tX = CliffordOperator(sX)
tY = CliffordOperator(sY)
tZ = CliffordOperator(sZ)
tH = tHadamard
tS = tPhase
tISWAP = (tI ⊗ tH) * tSWAP * tCNOT * tSWAP * tCNOT * (tH ⊗ tI) * (tS ⊗ tS)

stA = [tI, tH, tS, tH*tS, tS*tH, tH*tS*tH]
stB = [tI, tX, tY, tZ]
stR = [tI⊗tI, tSWAP, tISWAP, tCNOT]

"""Construct the `QuantumClifford` gate associated with a one-qubit Clifford index."""
function clifford_gate(ind::Clifford1Index)
    return stA[1+ind.a]*stB[1+ind.b]
end

"""Construct the `QuantumClifford` gate associated with a two-qubit Clifford index."""
function clifford_gate(ind::Clifford2Index)
    gate = clifford_gate(ind.c1) ⊗ clifford_gate(ind.c2)
    gate = stR[1+ind.r] * gate
    gate = (clifford_gate(ind.s1) ⊗ clifford_gate(ind.s2)) * gate
    return gate
end

function applyGate!(CPsi::Union{Stabilizer,CliffordOperator}, ind::Clifford1Index, n::Int)
    return applyGate!(CPsi, clifford_gate(ind), n)
end

function applyGate!(CPsi::Union{Stabilizer,CliffordOperator}, ind::Clifford2Index, n1::Int, n2::Int)
    return applyGate!(CPsi, clifford_gate(ind), n1, n2)
end

# -- Tensor Network Clifford Gates --------------------------------------------------------------------

tnA = ["", "H", "S", "HS", "SH", "HSH"]
tnB = ["", "X", "Y", "Z"]
tnR = ["I" "SWAP" "iSWAP" "CNOT"]

function clifford_gate(ind::Clifford1Index, index::Index)
    gates = ITensor[]
    gatenames = reverse(tnA[1+ind.a]*tnB[1+ind.b])
    for gatename in gatenames
        push!(gates,op(string(gatename),index))
    end
    operator = op("I", index)
    return ITensors.apply(gates, operator)
end

function clifford_gate(ind::Clifford2Index, index1::Index, index2::Index)
    gates = ITensor[]
    push!(gates,clifford_gate(ind.c1, index1))
    push!(gates,clifford_gate(ind.c2, index2))
    push!(gates,op(tnR[1+ind.r], index1, index2))
    push!(gates,clifford_gate(ind.s1, index1))
    push!(gates,clifford_gate(ind.s2, index2))
    operatorId = op("Id", index1) * op("Id", index2)
    return ITensors.apply(gates, operatorId)
end

function ITensors.op(ind::Clifford1Index, index::Index)
    return clifford_gate(ind, index)
end

function ITensors.op(ind::Clifford2Index, index1::Index, index2::Index)
    return clifford_gate(ind, index1, index2)
end

# -- Random Clifford Circuit ----------------------------------------------------------------------

"""Generate a random Clifford circuit on `N` qubits with the requested depth."""
function random_clifford_circuit(N::Int, debth::Int; clifford_set=Clifford2IndexSet(:all))
    @assert N>1
    res = one(CliffordOperator, N)
    for _ in 1:debth
        num = rand(1:length(clifford_set))-1
        ind = Clifford2Index(clifford_set, num)
        gate = clifford_gate(ind)
        sites = sort(StatsBase.sample(1:N, 2, replace=false))
        res = (gate, sites) * res
    end
    return res
end