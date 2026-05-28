using DisentangleCAMPS
using Test
using CliffordMPS

@testset "apply / ^ compatibility" begin
  for P in [P"XYZI", P"XXXX"]
    N = length(P)
    C = random_clifford_circuit(N, 2N^2)
    P_sum = PauliSum(Stabilizer(QuantumClifford.Tableau([P])))
    res_sum = P_sum^C
    @test res_sum.ops[1] == apply(P, C)
  end
end