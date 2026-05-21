using Revise

import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

# == Terminal input/output ==============================================================

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

# == File input/output ==================================================================

function save_three_columns(a, b, c, filename)
  n = max(length(a), length(b), length(c))
  open(filename, "a") do f
    for i in 1:n
      ai = i <= length(a) ? a[i] : ""
      bi = i <= length(b) ? b[i] : ""
      ci = i <= length(c) ? c[i] : ""
      println(f, "$(ai) $(bi) $(ci)")
    end
  end
end