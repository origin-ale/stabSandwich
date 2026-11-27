import PauliPropagation as pp
using QuantumClifford
using .DisentangleCAMPS
using Statistics

# == File input/output ==================================================================

""" ```stack_samples(sample_evs)```

Stack per-sample vectors into a layers × samples matrix, padding samples shorter than the \
longest one with ```missing``` at the end. The element type is promoted over all samples \
(e.g. an ```Int``` and a ```Float64``` sample give ```Float64```), so no value is \
converted to a narrower type."""
function stack_samples(sample_evs)
  maxlen = maximum(length, sample_evs)
  T = mapreduce(eltype, promote_type, sample_evs)
  evs = Matrix{Union{Missing, T}}(missing, maxlen, length(sample_evs))
  for (j, a) in enumerate(sample_evs)
    evs[1:length(a), j] = a
  end
  return evs
end

""" ```quantiles(xs)```

Lower quartile, median and upper quartile of ```xs``` (ordered for gnuplot \
```with candlesticks```). An empty ```xs``` gives NaNs."""
function quantiles(xs)
  isempty(xs) && return (NaN, NaN, NaN)
  q = quantile(xs, (0.25, 0.5, 0.75))
  return (q[1], q[2], q[3])
end

""" ```label_header(label, μ)```

Block header of the ```save_*``` writers: ```label = "n" => n``` gives \
```# n = <n>, μ = <μ>``` (count mode); the T11 replay passes ```"p" => p```, giving \
```# p = <p>, μ = <μ>```. With ```μ = nothing``` (models without a μ, e.g. \
```"ε" => ε``` for the kicked Ising model) the header is ```# <k> = <v>```. In the \
writers' docstrings, ```# n = …``` stands for any of these."""
label_header(label::Pair{<:AbstractString}, μ) = "# $(label.first) = $(label.second), μ = $μ"
label_header(label::Pair{<:AbstractString}, ::Nothing) = "# $(label.first) = $(label.second)"

""" ```save_full_samples(output_full, μ, label, evs)```

Append every sample's whole trajectory for one parameter point: a \
```# n = …, μ = …``` header, then for each column of ```evs``` (layers × samples, see \
```stack_samples```) a ```# sample s``` line and one ```layer<TAB>value``` line per \
non-```missing``` entry (layer = row − 1), followed by a blank line; the block ends with \
a second blank line (gnuplot block separator)."""
function save_full_samples(output_full, μ, label, evs)
  open(output_full, "a") do full_io
    println(full_io, label_header(label, μ))
    for sample_idx in axes(evs, 2)
      println(full_io, "# sample $sample_idx")
      for layer_idx in axes(evs, 1)
        ev = evs[layer_idx, sample_idx]
        ismissing(ev) && continue
        println(full_io, "$(layer_idx - 1)\t$(ev)")
      end
      println(full_io)
    end
    println(full_io)  # second blank line: gnuplot index (block) separator
  end
end

""" ```save_per_sample(output, μ, label, cols...)```

Append one block per parameter point holding per-sample scalars: a \
```# n = …, μ = …``` header, then one ```sample cols…``` row per sample, \
followed by two blank lines (gnuplot block separator)."""
function save_per_sample(output, μ, label, cols::AbstractVector...)
  open(output, "a") do io
    println(io, label_header(label, μ))
  end
  save_columns(output, collect(1:maximum(length, cols)), cols...)
  open(output, "a") do io
    print(io, "\n\n")  # two blank lines: gnuplot index (block) separator
  end
end

""" ```save_peaks(output, samples, μ, label)```

Append one ```sample peak``` row per sample, where the peak is the maximum \
over that sample's own evolution. Samples with no data (all ```missing```, \
e.g. a method that never reached its switch layer) are written as NaN. \
One block per parameter point, terminated by two blank lines."""
function save_peaks(output, samples, μ, label)
  isempty(samples) && return
  peaks = map(samples) do s
    # NaN sentinel rather than `init` so that integer resources stay integral
    any(!ismissing, s) ? maximum(skipmissing(s)) : NaN
  end
  save_per_sample(output, μ, label, peaks)
end

""" ```save_stats(output, evs, μ, label)```

Append the per-layer mean and standard error of ```evs``` (layers × samples, see \
```stack_samples```; ```missing``` entries skipped) as three space-separated columns \
(layer, mean, error), preceded by a ```# n = …, μ = …``` header line \
identifying the parameter point, and followed by two blank lines (gnuplot block \
separator). Layer = row − 1; rows with no data are left out. The error is the sample \
standard deviation (n − 1) over √n, NaN for a single sample."""
function save_stats(output, evs, μ, label)
  rows = [i for i in axes(evs, 1) if any(!ismissing, view(evs, i, :))]
  ev_means = [mean(skipmissing(evs[i, :])) for i in rows]
  ev_errs = [std(skipmissing(evs[i, :]))/sqrt(count(!ismissing, evs[i, :])) for i in rows]
  layers = rows .- 1

  open(output, "a") do io
    println(io, label_header(label, μ))
  end
  save_three_columns(layers, ev_means, ev_errs, output)
  open(output, "a") do io
    print(io, "\n\n")  # two blank lines: gnuplot index (block) separator
  end
end

""" ```save_stats_maxcol(output, samples, μ, label)```

Like ```save_stats```, but takes the raw per-sample resource vectors \
```samples``` and appends, as a fourth column sharing the same layer (cycle) \
column, the full per-layer evolution of the single sample that reached the \
highest resource value (peak over its own evolution), and as a fifth column \
the number of samples with data at each layer, and as five further columns \
the per-layer quartiles and extrema over the samples. The block header records \
the selected (1-based) sample index; columns are \
```layer, mean, error, max-sample, n-samples, q1, median, q3, min, max``` \
(columns 6-10 ordered for gnuplot ```with candlesticks```, whiskers at the \
extrema) and the block is followed by \
two blank lines (gnuplot block separator). Samples may contain leading \
```missing``` entries (e.g. resources of a method that starts mid-circuit): \
rows with no data are skipped and gaps in the max-sample column are written \
as NaN. Note that column 4 is one sample's whole trajectory, whereas column 10 \
is the per-layer maximum over all samples."""
function save_stats_maxcol(output, samples, μ, label)
  isempty(samples) && return
  evs = stack_samples(samples)
  rows = [i for i in axes(evs, 1) if any(!ismissing, view(evs, i, :))]
  counts = [count(!ismissing, view(evs, i, :)) for i in rows]
  ev_means = [mean(skipmissing(evs[i, :])) for i in rows]
  ev_errs = [std(skipmissing(evs[i, :]))/sqrt(count(!ismissing, evs[i, :])) for i in rows]
  quartiles = [quantiles(collect(skipmissing(view(evs, i, :)))) for i in rows]
  ev_q1 = [q[1] for q in quartiles]
  ev_medians = [q[2] for q in quartiles]
  ev_q3 = [q[3] for q in quartiles]
  ev_mins = [minimum(skipmissing(view(evs, i, :))) for i in rows]
  ev_maxs = [maximum(skipmissing(view(evs, i, :))) for i in rows]
  layers = rows .- 1

  imax = argmax([maximum(skipmissing(s); init = -Inf) for s in samples])
  maxcol = [coalesce(evs[i, imax], NaN) for i in rows]

  open(output, "a") do io
    println(io, label_header(label, μ), ", max sample = $imax")
  end
  save_columns(output, layers, ev_means, ev_errs, maxcol, counts,
    ev_q1, ev_medians, ev_q3, ev_mins, ev_maxs)
  open(output, "a") do io
    print(io, "\n\n")  # two blank lines: gnuplot index (block) separator
  end
end

""" ```save_doped_gates(output, gates, doped_inds_samples, layer_ends, μ, label)```

Append one block per parameter point listing the doped rotation gates of each \
sample: a ```# n = …, μ = …``` header, then for each sample a \
```# sample s: g doped gates``` line followed by one \
```gate-position cycle Pauli-string qubit-indices``` line per doped gate \
(e.g. ```17 2 XX 3 4```), where the position indexes the gates array and the \
(1-based) cycle is derived from ```layer_ends```. Samples are separated by a \
blank line and the block ends with two blank lines (gnuplot block separator)."""
function save_doped_gates(output, gates, doped_inds_samples, layer_ends, μ, label)
  open(output, "a") do io
    println(io, label_header(label, μ))
    for (s, inds) in enumerate(doped_inds_samples)
      println(io, "# sample $s: $(length(inds)) doped gates")
      for k in inds
        g = gates[k]
        cycle = searchsortedfirst(layer_ends, k)
        println(io, k, " ", cycle, " ",
          join(string.(g.symbols)), " ", join(g.qinds, " "))
      end
      println(io)
    end
    println(io)  # second blank line: gnuplot index (block) separator
  end
end

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

""" ```save_three_columns(a, b, c, filename)```

```save_columns(filename, a, b, c)```, with the file name last (the older argument order)."""
save_three_columns(a, b, c, filename) = save_columns(filename, a, b, c)