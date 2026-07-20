module CampsPP

include("CliffordMPSPauliPropagation.jl")
export paulivec
export getpauli
export leftover_rotgates

include("Utility.jl")
export stringtopauli_sym
export save_three_columns
export save_rows
export save_columns
export stack_samples
export save_full_samples
export save_per_sample
export save_peaks
export save_stats
export append_stats
export save_stats_maxcol
export save_doped_gates

include("Circuits.jl")
export dopeT
export dopeMagic
export sub_magic
export xy_magic
export z_magic
export rotation_circuit
export xxz_alt_circuit
export xxz_layer_circuit
export xxz_circuit
export fSim_circuit

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
