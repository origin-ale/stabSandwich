using Revise

using .DisentangleCAMPS
import PauliPropagation as pp
import MagicTensors as mt
using ITensors, ITensorMPS

using ProgressMeter

# == File manipulation ==================================================================

function initialize_output(output::AbstractString, obsname::AbstractString, params::AbstractDict)
  open(output, "w") do f
    pairs = ["$(string(k))=$(string(v))" for (k,v) in params]
    header = isempty(pairs) ? "" : (join(pairs, " ") * " ")
    println(f, "# $(header) \n# Cycle $obsname\n\n")
  end
end

function append_datapoint(output::AbstractString, x::Number, y::Number)
    open(output, "a") do f
        println(f, "$x\t$y")
    end
end

# == State types ========================================================================

nqubits_state(ψ::mt.QubitCAMPS) = mt.nsites(ψ)

"Convert a PauliPropagation `PauliSum` to the `MagicTensors.QubitPauliSum` that `expectation_value` accepts. Call once, before the time loop, in CAMPS evolution."
convert_observable(::mt.QubitCAMPS, obs::pp.PauliSum) = mt.QubitPauliSum(obs)

expectation_value(ψ::mt.QubitCAMPS, obs) = mt.expectation(ψ, obs)
expectation_value(ψ::ITensorMPS.MPS, obs::mt.QubitPauliSum) = mt.expectation(ψ, obs)

# == Append expectation values ==========================================================

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

# == ProgressMeter showvalues ===========================================================

showvalues_χ(χ, bd, N, k) = () -> [("Bond dimension (max $χ)", bd), ("Magic qubits (max $N)", k)]
showvalues_P(Pmax, P) = () -> [("Pauli terms (max $Pmax)", P)]
showvalues_ee(ee) = () -> [("Entanglement entropies", round.(ee; digits=2))]

# == Utility ============================================================================

function getlayer(i, layer_ends)
  isnothing(layer_ends) && return i
  for l in eachindex(layer_ends)
    layer_ends[l] >= i && return l
  end
end

# Layer-end disentangler of `campssrc_circuit_dynamics`: up to N sweeps of `strategy` with
# `criterion`, stopping once a sweep gains ≤ 6, with no SVD truncation. On these
# nearest-neighbour sweeps `:chi3` never applies a gate (see DisentangleQubitCAMPS.jl).
const LAYER_MIN_DIFF = 6
layer_disentangler(::mt.QubitCAMPS, criterion, strategy, N) =
  sweep_disentangler(criterion, strategy, N; min_diff = LAYER_MIN_DIFF)
layer_disentangle!(ψ::mt.QubitCAMPS, de) = mt.disentangle!(de, ψ)

meanweight(psum::pp.PauliSum) = mean(pp.countweight(psum))
meancoeff(psum::pp.PauliSum) = mean(abs, pp.coefficients(psum))

"Σ|c|² over the terms of `psum`."
sumsquares(psum::pp.PauliSum) = sum(abs2, pp.coefficients(psum); init = 0.0)

const PP_NORM_RTOL = 1e-8

"""
    check_pp_norm(norm_obs, paulisum; rtol = PP_NORM_RTOL)

Throw if Pauli propagation raised Σ|c|² above `norm_obs`, its value for the input
observable, by more than `rtol` relative. Pauli rotations and Cliffords preserve Σ|c|² and
truncation only removes terms, so a larger value means the propagated sum is corrupt (as
with the PauliPropagation 0.7.3 hash bug described at `__init__` in CampsPP.jl). Throwing
keeps its expectation values out of the output. Rounding drift is far below `rtol`.
"""
function check_pp_norm(norm_obs::Real, paulisum::pp.PauliSum; rtol::Real = PP_NORM_RTOL)
  norm_out = sumsquares(paulisum)
  norm_out > (1 + rtol) * norm_obs && error(
    "Pauli propagation raised Σ|c|² from $norm_obs to $norm_out (relative tolerance \
$rtol): the propagated sum is corrupt")
  return nothing
end

# == Circuit evolution ==================================================================

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
  progress = ProgressUnknown(desc = "Evolving CAMPS… gate ", enabled = showprogress)
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

    bd = DisentangleCAMPS.bonddim(ψ)
    next!(progress; showvalues = showvalues_χ(χ, bd, N, k))
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(bond dim. ≥ $χ)"
      println(f, "# CAMPS stopped at gate $i ", reason, "\n\n")
    end
  end

  track && return ψ, k, i, evs_camps, bonddims
  return ψ, k, i, evs_camps
end

"""```campssrc_circuit_dynamics(ψ, gates, phases, thl, obs[, output]; criterion = :entangle, strategy = :brickwork, showprogress = false, layer_ends = nothing, track = false)```

Evolve a copy of the CAMPS ψ through the whole circuit with the truncating `apply!(ψ, P, ϕ, thl)`.
⟨obs⟩ is recorded at the start and at every record point (the end of each layer in `layer_ends`,
or every gate if it is `nothing`), and also appended to the file `output` if one is given. At each
record point the layer-end disentangler (`criterion`, `strategy`, at most N sweeps, `min_diff = 6`)
runs until the run is flagged as truncated, i.e. once norm² < (1 − thl)·norm²(t = 0); after that it
no longer runs. Returns `(ψ, i, evs)`, plus the bond dimension at each record point if `track`.

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
track = false)
  return campssrc_circuit_dynamics(
    ψ_ext, gates, phases, thl, obs, nothing;
    criterion = criterion, strategy = strategy,
    showprogress = showprogress, layer_ends = layer_ends, track = track)
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
track = false)

  ψ = deepcopy(ψ_ext)
  N = nqubits_state(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs_camps = []
  bonddims = Int[]
  progress = ProgressUnknown(desc = "Evolving CAMPS (SVD)… gate ", enabled = showprogress)
  obs_mt = convert_observable(ψ, obs)
  layer_de = layer_disentangler(ψ, criterion, strategy, N)

  # Gates are unitary and disentangling Cliffords are exact, so in exact arithmetic the norm
  # only shrinks when the SVD cutoff in apply! discards weight. In floating point it also drifts
  # by rounding (about 1e-14), so for thl ≲ 1e-13 this check can fire with nothing truncated
  # (see the docstring warning).
  norm0_sq = norm(ψ)^2
  truncated = false

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
    end

    bd = DisentangleCAMPS.bonddim(ψ)
    next!(progress; showvalues = () -> [("Bond dimension", bd)])
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      println(f, "# CAMPS (SVD) stopped at gate $i (end of circuit)\n\n")
    end
  end

  track && return ψ, i, evs_camps, bonddims
  return ψ, i, evs_camps
end

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
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

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

    next!(progress; showvalues = showvalues_P(Nmax, NP))
  end

  finish!(progress)
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
  progress = ProgressUnknown(desc = "Evolving with Pauli prop… gate ", enabled = showprogress)

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

    next!(progress; showvalues = showvalues_P(Nmax, NP))
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = (i == M) ? "(end of circuit)" : "(n. of Paulis ≥ $Nmax)"
      println(f, "# Pauli prop. stopped at gate $i ", reason, "\n\n")
    end
  end
  track && return i, evs_pp, nterms, avgweights, avgcoeffs
  return i, evs_pp
end

function mps_circuit_dynamics(
  ψ_ext::MPS,
  gates::Vector,
  phases::Vector,
  thl,
  obs::pp.PauliSum,
  output;
  showprogress = false,
  layer_ends = nothing
)
  ψ = deepcopy(ψ_ext)                              
  N = length(ψ)
  M = length(gates)
  i = 0
  layer = 0
  evs = []
  progress = ProgressUnknown(desc = "Evolving MPS… gate ", enabled = showprogress)
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
    end
    # Entropies only for the progress display. mt.entanglement_entropy works on a copy, so
    # the gauge of ψ (and hence what the next truncation drops) does not depend on
    # showprogress. The old cmps.eEntropys! re-gauged ψ in place after every gate.
    if showprogress
      next!(progress; showvalues = showvalues_ee(mt.entanglement_entropy(ψ)))
    else
      next!(progress)
    end
  end

  finish!(progress)
  if !isnothing(output)
    open(output, "a") do f
      reason = "end of circuit"
      println(f, "# MPS stopped at gate $i ", reason, "\n\n")
    end
  end

  return ψ, evs
end