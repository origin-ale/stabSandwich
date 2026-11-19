using Test

include("helpers.jl")

@testset "CampsPP" begin
  include("TestCircuits.jl")
  include("TestInitialize.jl")
  include("TestTimeEvolution.jl")
  include("TestCliffordLargeN.jl")
  include("TestPPKeys.jl")
  include("TestMagicTensorsPP.jl")    # uses helpers from TestPPKeys.jl
  include("TestPerTermSeam.jl")       # post-port item 4; item 4 groups broken until the code lands
  include("TestPPHashBug.jl")
  include("TestDepwarn.jl")
  include("TestThreads.jl")
  include("TestUtilityFormat.jl")
  include("TestProgress.jl")          # uses helpers from TestTimeEvolution.jl
end
