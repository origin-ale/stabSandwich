using CampsPP.DisentangleCAMPS

using Test

# One wrapping testset: a failing file no longer aborts the run and hides the rest.
@testset "DisentangleCAMPS" begin
  # The self-test runs first, before QuantumClifford is loaded into this module: both
  # export `expect`, and the self-test uses the dense reference's unqualified.
  @testset "Dense reference (self-test)" begin
    include("TestDenseReference.jl")
  end
  include("helpers.jl")
  @testset "Golden set" begin
    include("TestGoldenSet.jl")
  end
  @testset "Disentangler" begin
    include("TestDisentangler.jl")
  end
  @testset "Evolve" begin
    include("TestEvolve.jl")
  end
  @testset "Dense checks" begin
    include("TestApplyDense.jl")
  end
  @testset "Clifford circuits at N = 46" begin
    include("TestCliffordLargeN.jl")
  end
  @testset "Entropies" begin
    include("TestEntropies.jl")
  end
  @testset "Layer-end disentangler (QubitCAMPS)" begin
    include("TestDisentangleQubitCAMPS.jl")
  end
end