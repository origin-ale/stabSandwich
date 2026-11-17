module DisentangleCAMPS

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
