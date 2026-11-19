using Test

include("helpers.jl")

@testset "CampsPP" begin
  include("TestCircuits.jl")
  include("TestInitialize.jl")
  include("TestTimeEvolution.jl")
  include("TestCliffordLargeN.jl")
  include("TestPPKeys.jl")
  include("TestMagicTensorsPP.jl")    # uses helpers from TestPPKeys.jl
  include("TestPPHashBug.jl")
  include("TestDepwarn.jl")
  include("TestThreads.jl")
  include("TestUtilityFormat.jl")
  include("TestProgress.jl")          # uses helpers from TestTimeEvolution.jl
end
