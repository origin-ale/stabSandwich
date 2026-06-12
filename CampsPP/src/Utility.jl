using Revise

import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS
using Statistics

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

function stack_samples(sample_evs)
  maxlen = maximum(length, sample_evs)
  T = eltype(first(sample_evs))
  evs = Matrix{Union{Missing, T}}(missing, maxlen, length(sample_evs))
  for (j, a) in enumerate(sample_evs)
    evs[1:length(a), j] = a
  end
  return evs
end

function save_full_samples(output_full, μ, evs)
  open(output_full, "a") do full_io
    println(full_io, "# μ = $μ")
    for sample_idx in axes(evs, 2)
      println(full_io, "# sample $sample_idx")
      for layer_idx in axes(evs, 1)
        ev = evs[layer_idx, sample_idx]
        ismissing(ev) && continue
        println(full_io, "$(layer_idx - 1)\t$(ev)")
      end
      println(full_io)
    end
  end
end

function save_stats(output, evs)
  ev_means = [mean(skipmissing(row)) for row in eachrow(evs)]
  ev_errs = [std(skipmissing(row))/sqrt(count(!ismissing, row)) for row in eachrow(evs)]
  layers = collect(0:size(evs, 1) - 1)

  save_three_columns(layers, ev_means, ev_errs, output)
  save_three_columns(["\n\n"], [""], [""], output)
end

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