"""
CAMPS evolution with disentangling on MagicTensors `QubitCAMPS` (submodule of CampsPP,
formerly a separate package): Pauli rotations with the Liu and Clark (2025) disentangler
(`apply!`, `evolve`, `evolve_deepcliffords`), the layer-end sweep disentangler
(`sweep_disentangler`) and plain-text progress lines (`ProgressLog`).
"""
module DisentangleCAMPS

# Plain-text progress lines (replaces ProgressMeter; used by the evolution loops here and in
# CampsPP)
include("Progress.jl")
export ProgressLog, isdue, needs_final, log_progress!, log_done!, log_message

include("Disentangler.jl")
export paulinature
export disentangler

# MagicTensors QubitCAMPS evolution (apply!, evolve, evolve_bonddim, evolve_deepcliffords)
include("EvolveQubitCAMPS.jl")
export evolve_deepcliffords, evolve
export rotation_gate

# Layer-end disentangler for QubitCAMPS: the old :chi3/:entangle × :snake/:brickwork with
# min_diff, on upstream disentangler types (plan, Phase 3)
include("DisentangleQubitCAMPS.jl")
export sweep_disentangler

end
