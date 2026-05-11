import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

function stringtopauli(str::String)
    char_to_sym = Dict('I' => :I, 'X' => :X, 'Y' => :Y, 'Z' => :Z)
    nq = length(str)
    paulis = Symbol[]
    qinds = Int[]
    for (i, c) in enumerate(str)
        s = char_to_sym[c]
        if s !== :I
            push!(paulis, s)
            push!(qinds, i)
        end
    end
    return pp.PauliString(nq, paulis, qinds)
end

""" ```random_rotations(t, N)```

Generate t N-qubit random Pauli rotation gates e^(iϕP), with phases ϕ returned separately.
"""
function random_rotations(t::Integer, Nqubits::Integer)
  paulistrings = [pp.inttosymbol(rand(0:4^Nqubits-1), Nqubits) for _ in 1:t]
  phases = 2π * rand(Float64, (t,)) # Exponential phases, ie. -1/2 * rotation angles
  qinds = [collect(1:Nqubits) for _ in paulistrings]
  rotations = pp.PauliRotation.(paulistrings, qinds)
  return rotations, phases
end

""" ```random_paulistring(N)```

Generate a random N-qubit Pauli string.
"""
random_paulistring(N::Integer) = stringtopauli(join(rand(["I","X","Y","Z"], N)))