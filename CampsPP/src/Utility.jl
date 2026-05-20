using Revise

import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

""" ```stringtopauli_sym(str::String)```

Convert a string consisting of ```I```, ```X```, ```Y```, ```Z``` \
to a vector of symbols ```:I```, ```:X```, ```:Y```, ```:Z``` \
(compatible with PauliPropagation.jl)."""
function stringtopauli_sym(str::String)
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

