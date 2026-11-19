"""
CAMPS-PP: Clifford-augmented MPS (MagicTensors `QubitCAMPS`) joined with Pauli propagation
(PauliPropagation.jl). Circuits are Pauli rotations exp(iϕP) given as `pp.PauliRotation`s
and phases ϕ (`xxz_circuit`, `rotation_circuit`, doped with `xy_magic`/`z_magic`); the
`*_circuit_dynamics` functions evolve a state (`computationalcamps`, `domainwallstate`) and
record an observable (e.g. `transferredmagnetization`) per layer (`layerends`), and the
`save_*` writers produce the gnuplot-ready output files. Conventions: test/CONVENTIONS.md.
The CAMPS evolution and disentangling live in the submodule `DisentangleCAMPS`.
"""
module CampsPP

# Formerly the separate package DisentangleCAMPS, merged in unchanged as a submodule.
include("DisentangleCAMPS.jl")
using .DisentangleCAMPS
export ProgressLog, isdue, needs_final, log_progress!, log_done!, log_message

# PauliPropagation ↔ MagicTensors bridge (converter, expectation seam, pp helpers)
include("MagicTensorsPauliPropagation.jl")
export getpauli
export leftover_rotgates
export pauliop_from_ppterm

include("Utility.jl")
export quantiles
export save_three_columns
export save_rows
export save_columns
export stack_samples
export save_full_samples
export save_per_sample
export save_peaks
export save_stats
export save_stats_maxcol
export save_doped_gates

include("Circuits.jl")
export doping_mask
export xy_magic
export z_magic
export rotation_circuit
export xxz_circuit

include("TimeEvolution.jl")
export camps_circuit_dynamics
export campssrc_circuit_dynamics
export pauliprop_circuit_dynamics
export campspp_circuit_dynamics
export mps_circuit_dynamics
export initialize_output
export append_datapoint

include("Initialize.jl")
export domainwallstate
export computationalcamps
export transferredmagnetization
export layerends

end
