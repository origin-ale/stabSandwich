using CampsPP
import CliffordMPS as cmps
import PauliPropagation as pp
using QuantumClifford
using DisentangleCAMPS

tCX = tCNOT
tCY = (tId1 ⊗ tPhase) * tCNOT * inv(tId1 ⊗ tPhase)
tCZ = (tId1 ⊗ tHadamard) * tCNOT * (tId1 ⊗ tHadamard)

k = 2
pauli = P"XIZYY"

expected_disentangler  = one(CliffordOperator, length(pauli))
apply!(expected_disentangler, tCX, [3,1])
apply!(expected_disentangler, tCZ, [3,4])
apply!(expected_disentangler, tCY, [3,5])
apply!(expected_disentangler, tSWAP, [3,4])
apply!(expected_disentangler, tPhase, [4])


# println("Creating disentangler for $pauli with $k magic qubits")
D, sign = disentangler(k, pauli)

println("=== Expected disentangler ===")
println(expected_disentangler)

println("=== Computed disentangler ===")
println(D)

same = (expected_disentangler == D) ? "coincide" : "DO NOT coincide"
println("\nThe two $(same).")

println("")