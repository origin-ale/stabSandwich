using .DisentangleCAMPS
import PauliPropagation as pp
import MagicTensors as mt
using ITensors, ITensorMPS


# == File manipulation ==================================================================

"""```initialize_output(output, obsname, params)```

Create (or overwrite) ```output``` with a two-line header, ```# k1=v1 k2=v2 …  ``` and \
```# Cycle <obsname>```, followed by two blank lines. ```params``` is a ```NamedTuple``` or a \
vector of ```key => value``` pairs, written in their own order, or an ```AbstractDict```, \
written in its iteration order (not stable across key sets or Julia versions). Note that a \
vector literal promotes its values (```[:N => 46, :t => 3.0]``` writes ```N=46.0```); pass \
strings or a ```NamedTuple``` to keep each value as given."""
function initialize_output(output::AbstractString, obsname::AbstractString, params)
  # Built before opening, so that unsupported `params` throw without truncating `output`
  kv = ["$(string(k))=$(string(v))" for (k, v) in header_pairs(params)]
  header = isempty(kv) ? "" : (join(kv, " ") * " ")
  open(output, "w") do f
    println(f, "# $(header) \n# Cycle $obsname\n\n")
  end
end

# The (key, value) pairs of `initialize_output`'s `params`; other types are a MethodError
header_pairs(params::NamedTuple) = pairs(params)
header_pairs(params::AbstractDict) = params
header_pairs(params::AbstractVector{<:Pair}) = params

"""```append_datapoint(output, x, y)```

Append the line ```x<TAB>y``` to the file ```output```. The ```*_circuit_dynamics``` \
functions write each recorded value this way, as ```record-point-index<TAB>value```."""
function append_datapoint(output::AbstractString, x::Number, y::Number)
    open(output, "a") do f
        println(f, "$x\t$y")
    end
end

# == State types ========================================================================

"Number of qubits of the state."
nqubits_state(ψ::mt.QubitCAMPS) = mt.nsites(ψ)

"Convert a PauliPropagation `PauliSum` to the `MagicTensors.QubitPauliSum` that `expectation_value` accepts. Call once, before the time loop, in CAMPS evolution."
convert_observable(::mt.QubitCAMPS, obs::pp.PauliSum) = mt.QubitPauliSum(obs)

"""```expectation_value(ψ, obs)```

⟨ψ|obs|ψ⟩, not normalised, possibly complex. On a `QubitCAMPS`, `obs` is a \
`MagicTensors.QubitPauliSum` (see `convert_observable`) or a PauliPropagation `PauliSum` \
(through the `mt.expectation` seam in MagicTensorsPauliPropagation.jl); on an MPS, a \
`QubitPauliSum`."""
expectation_value(ψ::mt.QubitCAMPS, obs) = mt.expectation(ψ, obs)
expectation_value(ψ::ITensorMPS.MPS, obs::mt.QubitPauliSum) = mt.expectation(ψ, obs)

# == Append expectation values ==========================================================

"""```append_expectation!(out_array, out_file, state, obs, i)```
```append_expectation!(out_array, out_file, ev, i)```

Record one value: push ```real(ev)``` to ```out_array``` and, unless ```out_file``` is \
```nothing```, append ```i<TAB>real(ev)``` to it. The value is \
```expectation_value(state, obs)``` for a `QubitCAMPS` or MPS ```state```, or the overlap \
of ```obs``` with the computational basis state that has 1s at ```state``` when \
```state``` is a vector of qubit indices (PauliPropagation's ```overlapwithcomputational```)."""
function append_expectation!(
out_array::Vector, out_file::AbstractString, ψ::mt.QubitCAMPS, obs, i::Integer)
  ev = expectation_value(ψ, obs)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, ::Nothing, ψ::mt.QubitCAMPS, obs, i::Integer)
  ev = expectation_value(ψ, obs)
  push!(out_array, real(ev))
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, ψ::ITensorMPS.MPS,
obs::mt.QubitPauliSum, i::Integer)
  ev = expectation_value(ψ, obs)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, ::Nothing, ψ::ITensorMPS.MPS,
obs::mt.QubitPauliSum, i::Integer)
  ev = expectation_value(ψ, obs)
  push!(out_array, real(ev))
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, onebitinds::Vector{<:Integer}, obs, i::Integer)
  ev = pp.overlapwithcomputational(obs, onebitinds)
  append_expectation!(out_array, out_file, ev, i)
end

function append_expectation!(
out_array::Vector, ::Nothing, onebitinds::Vector{<:Integer}, obs, i::Integer)
  ev = pp.overlapwithcomputational(obs, onebitinds)
  push!(out_array, real(ev))
end

function append_expectation!(
out_array::Vector, out_file::AbstractString, ev::Number, i::Integer)
  push!(out_array, real(ev))
  append_datapoint(out_file, i, real(ev))
end

function append_expectation!(
out_array::Vector, ::Nothing, ev::Number, i::Integer)
  push!(out_array, real(ev))
end

# == Progress values ====================================================================

# Values shown on the MPS progress lines: maximum bond dimension and maximum bond
# entanglement entropy (costs an SVD sweep, on a copy of ψ)
mps_progress_values(ψ::MPS) = ("bond dim." => maxlinkdim(ψ),
  "max S_E" => round(maximum(mt.entanglement_entropy(ψ); init = 0.0); digits = 2))

# == Utility ============================================================================

"""```getlayer(i, layer_ends)```

Index of the record point at or after gate ```i```: the first ```l``` with \
```layer_ends[l] ≥ i```, i.e. the layer that gate ```i``` belongs to, or ```i``` itself \
when ```layer_ends``` is ```nothing``` (every gate a record point). ```nothing``` if \
```i``` is past the last layer end."""
function getlayer(i, layer_ends)
  isnothing(layer_ends) && return i
  for l in eachindex(layer_ends)
    layer_ends[l] >= i && return l
  end
end

# Layer-end disentangler of `campssrc_circuit_dynamics`: up to N sweeps of `strategy` with
# `criterion`, stopping once a sweep gains ≤ 6, with no SVD keywords: the applied gates drop
# only numerical zeros (σ ≤ τ·σ_max, see DisentangleQubitCAMPS.jl and CONVENTIONS.md §2).
const LAYER_MIN_DIFF = 6
layer_disentangler(::mt.QubitCAMPS, criterion, strategy, N) =
  sweep_disentangler(criterion, strategy, N; min_diff = LAYER_MIN_DIFF)
layer_disentangle!(ψ::mt.QubitCAMPS, de) = mt.disentangle!(de, ψ)

"Mean Pauli weight (number of non-identity sites) over the terms of `psum`."
meanweight(psum::pp.PauliSum) = mean(pp.countweight(psum))
"Mean |c| over the terms of `psum`."
meancoeff(psum::pp.PauliSum) = mean(abs, pp.coefficients(psum))

"Σ|c|² over the terms of `psum`."
sumsquares(psum::pp.PauliSum) = sum(abs2, pp.coefficients(psum); init = 0.0)

const PP_NORM_RTOL = 1e-8

"""
    check_pp_norm(norm_obs, paulisum; rtol = PP_NORM_RTOL)

Throw if Pauli propagation raised Σ|c|² above `norm_obs`, its value for the input
observable, by more than `rtol` relative. Pauli rotations and Cliffords preserve Σ|c|² and
truncation only removes terms, so a larger value means the propagated sum is corrupt (as
with the PauliPropagation 0.7.3 hash bug, fixed by PauliPropagation 0.9: 0.7.3 created
Pauli-string integer types at runtime, and a `PauliSum` filled in the world age where its
type was created was corrupted once used in a later one). Throwing keeps its expectation
values out of the output. Rounding drift is far below `rtol`.
"""
function check_pp_norm(norm_obs::Real, paulisum::pp.PauliSum; rtol::Real = PP_NORM_RTOL)
  norm_out = sumsquares(paulisum)
  norm_out > (1 + rtol) * norm_obs && error(
    "Pauli propagation raised Σ|c|² from $norm_obs to $norm_out (relative tolerance \
$rtol): the propagated sum is corrupt")
  return nothing
end

# == Circuit evolution ==================================================================
#
# Common to all `*_circuit_dynamics` functions (test/CONVENTIONS.md §2):
# - `gates` are `pp.PauliRotation`s and `phases` the ϕ of each gate exp(iϕP) (PP itself is
#   given the angle −2ϕ, as `pp.PauliRotation(θ)` is exp(−iθP/2));
# - ⟨obs⟩ is recorded at t = 0 and at every record point: the gate indices in `layer_ends`,
#   or every gate if `layer_ends = nothing`. Entry l + 1 of the result is the value after
#   record point l. Values are real parts, not normalised by ⟨ψ|ψ⟩;
# - `output`, if given (not `nothing`), is a file to append each value to as
#   `l<TAB>value` (create it first with `initialize_output`), then a
#   `# … stopped at gate i (reason)` line and two blank lines;
# - the input state is copied, not modified;
# - `showprogress = true` prints plain-text progress lines (see `ProgressLog`).

"""```campspp_circuit_dynamics(ψ, χ, thl, Nmax, gates, phases, obs[, output]; showprogress = false, k = 0, layer_ends = nothing, track = false)```

CAMPS-PP ("sandwich") evolution of a copy of the CAMPS ψ: ```camps_circuit_dynamics``` \
until the bond dimension reaches ```χ```, then ```pauliprop_circuit_dynamics``` for the \
remaining gates, evaluating the propagated observable on the CAMPS state at the switch. \
PP drops terms with |c| < ```thl``` and stops after the first record point at which the \
propagated sum has ≥ ```Nmax``` terms. ```k``` is the number of magic qubits ψ already \
has. Recording, ```gates```/```phases``` and ```output``` as described at the top of this \
section (in TimeEvolution.jl); the CAMPS and PP parts write to the same file.

Returns ```(evs, s)```: one value per record point reached, starting at t = 0 (CAMPS \
values then PP values, with no gap or repeat at the switch), and ```s```, the last gate \
applied (```length(gates)``` unless PP stopped on ```Nmax```). With ```track = true```, \
returns ```(evs, s, bonddims, nterms, avgweights, avgcoeffs, s_switch, k_switch)```: the \
CAMPS bond dimension at t = 0 and at each CAMPS record point; the number of Pauli terms, \
mean weight and mean |c| of the observable at the switch and then at each PP record point; \
the last gate applied in CAMPS (PP starts at ```s_switch + 1```) and the number of magic \
qubits there."""
function campspp_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
χ::Integer,
thl::Real,
Nmax::Integer,
gates::Vector{},
phases::Vector,
obs::pp.PauliSum,
output;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)
  ψ = deepcopy(ψ_ext)

  if track
    # Keep the switch point: `s_switch` is the last gate applied in CAMPS
    # (PP starts at `s_switch + 1`; it equals length(gates) if CAMPS never hit
    # χ) and `k_switch` the number of magic qubits there
    ψ_evo, k_switch, s_switch, evs_camps, bonddims = camps_circuit_dynamics(
      ψ, gates, phases, χ, obs, output;
      showprogress = showprogress, k = k, layer_ends = layer_ends, track = true)

    s, evs_pp, nterms, avgweights, avgcoeffs = pauliprop_circuit_dynamics(
      ψ_evo, s_switch, gates, phases, thl, Nmax, obs, output;
      showprogress = showprogress, layer_ends = layer_ends, track = true)
  else
    ψ_evo, _, s, evs_camps = camps_circuit_dynamics(
      ψ, gates, phases, χ, obs, output;
      showprogress = showprogress, k = k, layer_ends = layer_ends)

    s, evs_pp = pauliprop_circuit_dynamics(
      ψ_evo, s, gates, phases, thl, Nmax, obs, output;
      showprogress=showprogress, layer_ends = layer_ends)
  end

  evs_tot = []
  append!(evs_tot, evs_camps)
  append!(evs_tot, evs_pp)

  track && return evs_tot, s, bonddims, nterms, avgweights, avgcoeffs,
    s_switch, k_switch
  return evs_tot, s
end

function campspp_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
χ::Integer,
thl::Real,
Nmax::Integer,
gates::Vector{},
phases::Vector,
obs::pp.PauliSum;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)
  return campspp_circuit_dynamics(
    ψ_ext, χ, thl, Nmax, gates, phases, obs, nothing;
    showprogress = showprogress, k = k, layer_ends = layer_ends, track = track)
end

"""```camps_circuit_dynamics(ψ, gates, phases, χ, obs[, output]; showprogress = false, k = 0, layer_ends = nothing, track = false)```

Evolve a copy of the CAMPS ψ, which has ```k``` magic qubits, gate by gate with \
```apply!(ψ, k, P, ϕ)``` (exact Cliffords, disentangling, no truncation), recording ⟨obs⟩ \
as described at the top of this section (in TimeEvolution.jl). Stops at the end of the \
circuit, or after the gate that brings the bond dimension to ≥ ```χ``` (checked before \
each gate). Values are recorded only at record points, so a stop inside a layer leaves \
that layer unrecorded.

Returns ```(ψ, k, i, evs)```: the evolved CAMPS, its number of magic qubits, the last gate \
applied and the recorded values; with ```track = true``` also the bond dimension at t = 0 \
and at each record point."""
function camps_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
gates::Vector,
phases::Vector,
χ::Integer,
obs::pp.PauliSum;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)
  return camps_circuit_dynamics(
    ψ_ext, gates, phases, χ, obs, nothing;
    showprogress = showprogress, k = k, layer_ends = layer_ends, track = track)
end

function camps_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
gates::Vector,
phases::Vector,
χ::Integer,
obs::pp.PauliSum,
output;
showprogress = false,
k = 0,
layer_ends = nothing,
track = false)

  ψ = deepcopy(ψ_ext)
  N = nqubits_state(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs_camps = []
  bonddims = Int[]
  progress = ProgressLog("CAMPS gate", M; enabled = showprogress, indent = 2)
  obs_mt = convert_observable(ψ, obs)

  append_expectation!(evs_camps, output, ψ, obs_mt, 0)
  track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
  layer += 1

  while DisentangleCAMPS.bonddim(ψ) < χ && i < M
    i += 1
    gate = PauliOperator(getpauli(gates[i], N))
    phase = phases[i]

    k = apply!(ψ, k, gate, phase)

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs_camps, output, ψ, obs_mt, layer)
      track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
      layer += 1
    end

    showprogress && log_progress!(progress, i,
      "bond dim." => (DisentangleCAMPS.bonddim(ψ), χ), "k" => (k, N))
  end

  showprogress && log_done!(progress, i,
    "bond dim." => (DisentangleCAMPS.bonddim(ψ), χ), "k" => (k, N))
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(bond dim. ≥ $χ)"
      println(f, "# CAMPS stopped at gate $i ", reason, "\n\n")
    end
  end

  track && return ψ, k, i, evs_camps, bonddims
  return ψ, k, i, evs_camps
end

"""```campssrc_circuit_dynamics(ψ, gates, phases, thl, obs[, output]; criterion = :entangle, strategy = :brickwork, showprogress = false, layer_ends = nothing, track = false, χmax = nothing)```

CAMPS with SVD truncation instead of disentangling ("CAMPS (SVD)" in the progress lines).
Evolve a copy of the CAMPS ψ through the whole circuit with the truncating `apply!(ψ, P, ϕ, thl)`.
⟨obs⟩ is recorded at the start and at every record point (the end of each layer in `layer_ends`,
or every gate if it is `nothing`), and also appended to the file `output` if one is given. At each
record point the layer-end disentangler (`criterion`, `strategy`, at most N sweeps, `min_diff = 6`)
runs until the run is flagged as truncated, i.e. once norm² < (1 − thl)·norm²(t = 0); after that it
no longer runs. With an integer `χmax`, the run stops after the first record point at which the bond
dimension (measured after the layer-end disentangler) is ≥ `χmax`; that value is still recorded. With
the default `nothing` it runs to the end of the circuit. Returns `(ψ, i, evs)`, plus the bond
dimension at each record point if `track`; `i` is the last gate applied.

!!! warning
    For thl ≲ 1e-13 the truncation flag is unreliable: the rounding drift of the norm (about
    1e-14) can set it when nothing has been truncated. Whether and when the layer-end disentangler
    runs then depends on floating-point details (it differed between the old CliffordMPS stack and
    MagicTensors). Use thl ≥ 1e-11 when the disentangler's behaviour matters."""
function campssrc_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
gates::Vector,
phases::Vector,
thl::Real,
obs::pp.PauliSum;
criterion = :entangle,
strategy = :brickwork,
showprogress = false,
layer_ends = nothing,
track = false,
χmax = nothing)
  return campssrc_circuit_dynamics(
    ψ_ext, gates, phases, thl, obs, nothing;
    criterion = criterion, strategy = strategy,
    showprogress = showprogress, layer_ends = layer_ends, track = track, χmax = χmax)
end

function campssrc_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
gates::Vector,
phases::Vector,
thl::Real,
obs::pp.PauliSum,
output;
criterion = :entangle,
strategy = :brickwork,
showprogress = false,
layer_ends = nothing,
track = false,
χmax = nothing)

  ψ = deepcopy(ψ_ext)
  N = nqubits_state(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs_camps = []
  bonddims = Int[]
  progress = ProgressLog("CAMPS (SVD) gate", M; enabled = showprogress, indent = 2)
  obs_mt = convert_observable(ψ, obs)
  layer_de = layer_disentangler(ψ, criterion, strategy, N)

  # Gates are unitary and disentangling Cliffords are exact, so in exact arithmetic the norm
  # only shrinks when the SVD cutoff in apply! discards weight. In floating point it also drifts
  # by rounding (about 1e-14), so for thl ≲ 1e-13 this check can fire with nothing truncated
  # (see the docstring warning).
  norm0_sq = norm(ψ)^2
  truncated = false
  stopped = false

  append_expectation!(evs_camps, output, ψ, obs_mt, 0)
  track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
  layer += 1

  while i < M
    i += 1
    gate = PauliOperator(getpauli(gates[i], N))
    phase = phases[i]

    apply!(ψ, gate, phase, thl)

    if !truncated && norm(ψ)^2 < norm0_sq * (1 - thl)
      truncated = true
    end

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs_camps, output, ψ, obs_mt, layer)
      truncated || layer_disentangle!(ψ, layer_de)
      track && push!(bonddims, DisentangleCAMPS.bonddim(ψ))
      layer += 1
      # Budget stop: only at record points, after the disentangler
      if !isnothing(χmax) && DisentangleCAMPS.bonddim(ψ) >= χmax
        stopped = true
      end
    end

    showprogress && log_progress!(progress, i, "bond dim." => DisentangleCAMPS.bonddim(ψ))
    stopped && break
  end

  showprogress && log_done!(progress, i, "bond dim." => DisentangleCAMPS.bonddim(ψ))
  if !isnothing(output)
    open(output, "a") do f
      reason = stopped ? "(bond dim. ≥ $χmax)" : "(end of circuit)"
      println(f, "# CAMPS (SVD) stopped at gate $i ", reason, "\n\n")
    end
  end

  track && return ψ, i, evs_camps, bonddims
  return ψ, i, evs_camps
end

"""```pauliprop_circuit_dynamics(ψ, start_gate, gates, phases, thl, Nmax, obs[, output]; showprogress = false, layer_ends = nothing, track = false)```
```pauliprop_circuit_dynamics(onebitinds, gates, phases, thl, Nmax, obs[, output]; showprogress = false, layer_ends = nothing, track = false)```

Pauli propagation (Heisenberg picture) of ```obs``` through the circuit, recording ⟨obs⟩ \
as described at the top of this section (in TimeEvolution.jl). At each record point the \
observable is propagated from scratch back through all PP gates so far (```pp.propagate``` \
with ```min_abs_coeff = thl```), checked with ```check_pp_norm```, and evaluated on the \
initial state. Stops at the end of the circuit or after the first record point at which \
the propagated sum has ≥ ```Nmax``` terms (the value there is still recorded).

- With a CAMPS ```ψ``` (the continuation of CAMPS-PP): the gates after ```start_gate``` are \
  propagated and evaluated on ψ (which is not evolved). Only the record points after \
  ```start_gate``` are returned, with no t = 0 entry.
- With ```onebitinds```: the initial state is the computational basis state with 1s at \
  those qubits, and the result starts with the t = 0 value.

Returns ```(i, evs)```, the last gate propagated and the recorded values; with \
```track = true```, ```(i, evs, nterms, avgweights, avgcoeffs)```: number of terms, mean \
Pauli weight and mean |c| of the unpropagated ```obs```, then of the propagated sum at each \
record point."""
function pauliprop_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
start_gate::Integer,
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum;
showprogress = false,
layer_ends = nothing,
track = false)
  return pauliprop_circuit_dynamics(
    ψ_ext, start_gate, gates, phases, thl, Nmax, obs, nothing;
    showprogress = showprogress, layer_ends = layer_ends, track = track)
end

function pauliprop_circuit_dynamics(
ψ_ext::mt.QubitCAMPS,
start_gate::Integer,
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum,
output;
showprogress = false,
layer_ends = nothing,
track = false)

  ψ = deepcopy(ψ_ext)
  M = length(gates)
  i = start_gate
  NP = 1
  layer = getlayer(start_gate+1, layer_ends)
  gates_pp = []
  angles_pp = []
  evs_pp = []
  nterms = Int[]
  avgweights = Float64[]
  avgcoeffs = Float64[]
  norm_obs = sumsquares(obs)
  progress = ProgressLog("PP gate", M; enabled = showprogress, indent = 2)

  track && push!(nterms, length(obs))
  track && push!(avgweights, meanweight(obs))
  track && push!(avgcoeffs, meancoeff(obs))

  while NP < Nmax && i < M
    i += 1
    gate = gates[i]
    angle = -2. * phases[i]

    push!(gates_pp, gate)
    push!(angles_pp, angle)

    if isnothing(layer_ends) || i == layer_ends[layer]
      paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
      check_pp_norm(norm_obs, paulisum)
      NP = length(paulisum)
      track && push!(nterms, NP)
      track && push!(avgweights, meanweight(paulisum))
      track && push!(avgcoeffs, meancoeff(paulisum))
      append_expectation!(evs_pp, output, ψ, paulisum, layer)
      layer += 1
    end

    showprogress && log_progress!(progress, i, "n. of Paulis" => (NP, Nmax))
  end

  showprogress && log_done!(progress, i, "n. of Paulis" => (NP, Nmax))
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
      println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
    end
  end
  track && return i, evs_pp, nterms, avgweights, avgcoeffs
  return i, evs_pp
end

function pauliprop_circuit_dynamics(
onebitinds::Vector{<:Integer},
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum;
showprogress = false,
layer_ends = nothing,
track = false)
  return pauliprop_circuit_dynamics(
    onebitinds, gates, phases, thl, Nmax, obs, nothing;
    showprogress = showprogress, layer_ends = layer_ends, track = track)
end

function pauliprop_circuit_dynamics(
onebitinds::Vector{<:Integer},
gates::Vector,
phases::Vector,
thl::Real,
Nmax::Integer,
obs::pp.PauliSum,
output::Union{AbstractString, Nothing};
showprogress = false,
layer_ends = nothing,
track = false)

  M = length(gates)
  i = 0
  NP = 1
  layer = 0
  gates_pp = []
  angles_pp = []
  evs_pp = []
  nterms = Int[]
  avgweights = Float64[]
  avgcoeffs = Float64[]
  norm_obs = sumsquares(obs)
  progress = ProgressLog("PP gate", M; enabled = showprogress, indent = 2)

  append_expectation!(evs_pp, output, onebitinds, obs, 0)
  track && push!(nterms, length(obs))
  track && push!(avgweights, meanweight(obs))
  track && push!(avgcoeffs, meancoeff(obs))
  layer += 1

  while NP < Nmax && i < M
    i += 1
    gate = gates[i]
    angle = -2. * phases[i]

    push!(gates_pp, gate)
    push!(angles_pp, angle)

    if isnothing(layer_ends) || i == layer_ends[layer]
      paulisum = pp.propagate(gates_pp, obs, angles_pp; min_abs_coeff = thl)
      check_pp_norm(norm_obs, paulisum)
      NP = length(paulisum)
      track && push!(nterms, NP)
      track && push!(avgweights, meanweight(paulisum))
      track && push!(avgcoeffs, meancoeff(paulisum))
      append_expectation!(evs_pp, output, onebitinds, paulisum, layer)
      layer += 1
    end

    showprogress && log_progress!(progress, i, "n. of Paulis" => (NP, Nmax))
  end

  showprogress && log_done!(progress, i, "n. of Paulis" => (NP, Nmax))
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
      println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
    end
  end
  track && return i, evs_pp, nterms, avgweights, avgcoeffs
  return i, evs_pp
end

"""```mps_circuit_dynamics(ψ, gates, phases, thl, obs, output; showprogress = false, layer_ends = nothing, χmax = nothing)```

Plain MPS evolution of a copy of the ITensors MPS ψ: each gate exp(iϕP) is applied as an \
MPO (```MagicTensors.apply!```, ITensorMPS densitymatrix algorithm) with SVD cutoff \
```thl```, and the state is renormalised after every gate. Records ⟨obs⟩ as described at \
the top of this section (in TimeEvolution.jl); ```output``` is required, pass \
```nothing``` for no file. With the default ```χmax = nothing``` it runs to the end of the \
circuit; with an integer ```χmax``` it stops after the first record point at which \
```maxlinkdim(ψ)``` is ≥ ```χmax``` (that value is still recorded). Returns \
```(ψ, evs)```."""
function mps_circuit_dynamics(
  ψ_ext::MPS,
  gates::Vector,
  phases::Vector,
  thl,
  obs::pp.PauliSum,
  output;
  showprogress = false,
  layer_ends = nothing,
  χmax = nothing
)
  ψ = deepcopy(ψ_ext)                              
  N = length(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs = []
  stopped = false
  progress = ProgressLog("MPS gate", M; enabled = showprogress, indent = 2)
  obs_mt = mt.QubitPauliSum(obs)

  append_expectation!(evs, output, ψ, obs_mt, 0)
  layer += 1

  while i < M
    i += 1
    gate = rotation_gate(PauliOperator(getpauli(gates[i], N)), phases[i]) # exp(iϕP)

    mt.apply!(ψ, gate; cutoff = thl) # MPO apply, ITensorMPS densitymatrix algorithm
    normalize!(ψ)

    if isnothing(layer_ends) || i == layer_ends[layer] # Works because || short circuits
      append_expectation!(evs, output, ψ, obs_mt, layer)
      layer += 1
      # Budget stop: only at record points
      if !isnothing(χmax) && maxlinkdim(ψ) >= χmax
        stopped = true
      end
    end
    # Entropies only for the progress display, and only when a line is due.
    # mt.entanglement_entropy works on a copy, so the gauge of ψ (and hence what the next
    # truncation drops) does not depend on showprogress.
    if showprogress && isdue(progress, i)
      log_progress!(progress, i, mps_progress_values(ψ)...)
    end
    stopped && break
  end

  # Entropies cost an SVD sweep: only when the final line is not printed yet
  if showprogress && needs_final(progress, i)
    log_done!(progress, i, mps_progress_values(ψ)...)
  end
  if !isnothing(output)
    open(output, "a") do f
      reason = stopped ? "(bond dim. ≥ $χmax)" : "(end of circuit)"
      println(f, "# MPS stopped at gate $i ", reason, "\n\n")
    end
  end

  return ψ, evs
end