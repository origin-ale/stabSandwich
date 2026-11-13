module DisentangleCAMPS

include("Disentangler.jl")
export paulinature
export disentangler

include("Evolve.jl")
export evolve_deepcliffords, evolve

# MagicTensors QubitCAMPS methods, staged next to the old CAMPS ones (plan, Phase 2)
include("EvolveQubitCAMPS.jl")
export rotation_gate

# Layer-end disentangler for QubitCAMPS: the old :chi3/:entangle × :snake/:brickwork with
# min_diff, on upstream disentangler types (plan, Phase 3)
include("DisentangleQubitCAMPS.jl")
export sweep_disentangler

end
