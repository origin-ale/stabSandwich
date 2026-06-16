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

function save_full_samples(output_full, μ, magic_prob, evs)
  open(output_full, "a") do full_io
    println(full_io, "# p = $magic_prob, μ = $μ")
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

""" ```save_stats(output, evs, μ, magic_prob)```

Append the per-layer mean and standard error of ```evs``` as three columns \
(layer, mean, error), preceded by a ```# magic_prob = …, μ = …``` header line \
identifying the parameter point, and followed by two blank lines (gnuplot block \
separator)."""
function save_stats(output, evs, μ, magic_prob)
  ev_means = [mean(skipmissing(row)) for row in eachrow(evs)]
  ev_errs = [std(skipmissing(row))/sqrt(count(!ismissing, row)) for row in eachrow(evs)]
  layers = collect(0:size(evs, 1) - 1)

  open(output, "a") do io
    println(io, "# p = $magic_prob, μ = $μ")
  end
  save_three_columns(layers, ev_means, ev_errs, output)
  open(output, "a") do io
    print(io, "\n\n")  # two blank lines: gnuplot index (block) separator
  end
end

""" ```append_stats(output, evs, μ, magic_prob)```

Alias for ```save_stats```, kept for backwards compatibility."""
append_stats(output, evs, μ, magic_prob) = save_stats(output, evs, μ, magic_prob)

""" ```save_rows(filename, params, rows; [blockend])```

Append one space-separated row per parameter point: \
```params[i] rows[i][1] rows[i][2] …```. If ```blockend = true```, terminate \
the block with a blank line (gnuplot-style block separator).
Always appends; call ```initialize_output``` first to (re)create the file."""
function save_rows(filename::AbstractString, params::AbstractVector,
                   rows::AbstractVector; blockend::Bool = false)
  open(filename, "a") do f
    for i in eachindex(params)
      println(f, join((params[i], rows[i]...), " "))
    end
    blockend && println(f)
  end
end

""" ```save_columns(filename, cols...)```

Append the given columns side by side, space-separated; shorter columns \
are padded with empty fields.
Always appends; call ```initialize_output``` first to (re)create the file."""
function save_columns(filename::AbstractString, cols::AbstractVector...)
  n = maximum(length, cols)
  open(filename, "a") do f
    for i in 1:n
      println(f, join((i <= length(c) ? string(c[i]) : "" for c in cols), " "))
    end
  end
end

save_three_columns(a, b, c, filename) = save_columns(filename, a, b, c)