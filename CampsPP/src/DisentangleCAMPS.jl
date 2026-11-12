module DisentangleCAMPS

include("Disentangler.jl")
export paulinature
export disentangler

include("Evolve.jl")
export evolve_deepcliffords, evolve

# MagicTensors QubitCAMPS methods, staged next to the old CAMPS ones (plan, Phase 2)
include("EvolveQubitCAMPS.jl")
export rotation_gate

end
