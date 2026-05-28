using Test

@testset "VectorizedMPDO" begin
    # even-length MPS should construct
    mps_even = random_qubit_mps(4, type="ITensors")
    @test VectorizedMPDO(mps_even) isa VectorizedMPDO

    # odd-length MPS should throw
    mps_odd = random_qubit_mps(3, type="ITensors")
    @test_throws ArgumentError VectorizedMPDO(mps_odd)
end
