# Phase 5 benchmarks: helpers shared by gen_sums.jl, bench_old.jl and bench_new.jl.
# Uses only functions that exist on both stacks (old `pre-port` worktree and `main`):
# `xxz_circuit`, `transferredmagnetization` and PauliPropagation, which the port left unchanged.

import PauliPropagation as pp
using Printf

flushprintln(args...) = (println(args...); flush(stdout))

const T8DIR = joinpath(@__DIR__, "..", "..", "CampsPP", "test", "data", "t8_handoff")
const SUMDIR = joinpath(@__DIR__, "..", "..", "output", "phase5", "sums")
const N46 = 46
const MAGIC_PHASE = 3π/16
const THL = 1e-10
const CHI = 64

"""The T8 hand-off sample (p = 0.06, sample 1): N, switch gate s, excited qubits, doped gates."""
function read_state(dir = T8DIR)
  fields = Dict(first(split(l, '\t')) => last(split(l, '\t'; limit = 2))
                for l in eachline(joinpath(dir, "state.txt")) if !startswith(l, '#'))
  N = parse(Int, fields["N"])
  s = parse(Int, fields["s"])
  excited = parse.(Int, split(fields["excited"]))
  doped = parse.(Int, split(fields["doped_gates"]))
  return (; N, s, excited, doped)
end

"""The doped XXZ circuit of the T8 sample (`xxz_circuit(π/4, π/4, 23, 46)`, doped angles 3π/16)."""
function sample_circuit(st, xxz_circuit)
  gates, phases = xxz_circuit(π/4, π/4, 23, st.N)
  phases[st.doped] .= MAGIC_PHASE
  return gates, phases
end

# Spec bit rule (test/CONVENTIONS.md §3): 2 bits per qubit, qubit j at bits 2(j−1), 2(j−1)+1.
termstring(term, N) = String(["IXYZ"[Int((term >> (2 * (j - 1))) & 0x3) + 1] for j in 1:N])

function ppsum_from_rows(N, rows)
  S = pp.PauliSum(N)
  for (str, c) in rows
    qs = [j for j in 1:N if str[j] != 'I']
    pp.add!(S, pp.PauliString(N, [Symbol(str[j]) for j in qs], qs, c))
  end
  return S
end

"""Read a sum file (T8 format: `# key\\tvalue` header lines, then `string\\tcoeff` rows)."""
function read_sum(path)
  header = Dict{String, String}()
  rows = Tuple{String, Float64}[]
  for l in eachline(path)
    if startswith(l, '#')
      startswith(l, "# ") && occursin('\t', l) &&
        (kv = split(l[3:end], '\t'; limit = 2); header[kv[1]] = kv[2])
    else
      str, c = split(l, '\t')
      push!(rows, (String(str), parse(Float64, c)))
    end
  end
  N = length(first(rows)[1])
  S = ppsum_from_rows(N, rows)
  haskey(header, "nterms") && @assert length(S) == parse(Int, header["nterms"])
  return S, header, rows
end

"""The benchmark sums, smallest first: the two T8 sums plus the generated ones in SUMDIR."""
function sum_files()
  files = [joinpath(T8DIR, "sum_0201.txt"), joinpath(T8DIR, "sum_0538.txt")]
  isdir(SUMDIR) && append!(files, sort(filter(f -> startswith(basename(f), "sum_"),
                                              readdir(SUMDIR; join = true))))
  nterms(f) = parse(Int, match(r"sum_0*(\d+)\.txt", basename(f))[1])
  return sort(files; by = nterms)
end
