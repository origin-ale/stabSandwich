module CampsPP

include("CiffordMPSPauliPropagation.jl")
export paulivec
export getpauli
export leftover_rotgates

include("Utility.jl")
export stringtopauli_sym
export save_three_columns

include("Circuits.jl")
export dopeT
export dopeMagic
export subMagic
export rotation_circuit
export xxz_layer_circuit
export xxz_circuit
export fSim_circuit

include("TimeEvolution.jl")
export camps_circuit_dynamics
export pauliprop_circuit_dynamics
export campspp_circuit_dynamics
export initialize_output

include("Initialize.jl")
export domainwallstate
export computationalcamps
export transferredmagnetization
export layerends

end
