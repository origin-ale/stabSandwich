using Revise

using CampsPP
import PauliPropagation as pp
import CliffordMPS as cmps

Nstrs = 10
Nq = 7 # Can handle up to 16384 terms
# Nq = Int(ceil(log(4, Nstrs+1)))
print("Conversion of $Nstrs strings on $Nq qubits ")

testdict = Dict(Int(i) => 1/Nstrs for i = 1:Nstrs)
testsum = pp.PauliSum(Nq, testdict)

_, ct, __ = @timed cmps.PauliSum(testsum)
print("took $ct s.")