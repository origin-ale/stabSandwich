using Test

include("helpers.jl")

@testset "CampsPP" begin
  include("TestCircuits.jl")
  include("TestInitialize.jl")
  include("TestTimeEvolution.jl")
  include("TestCliffordLargeN.jl")
  include("TestMagicTensorsPP.jl")
  include("TestPPHashBug.jl")
  include("TestUtilityFormat.jl")
end
