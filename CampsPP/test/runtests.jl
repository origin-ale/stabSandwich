using Test

include("helpers.jl")

@testset "CampsPP" begin
  include("TestCliffordMPSPauliPropagation.jl")
  include("TestCircuits.jl")
  include("TestInitialize.jl")
  include("TestTimeEvolution.jl")
end
