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
