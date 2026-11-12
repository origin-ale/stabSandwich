# Test entry point of CampsPP (run with `Pkg.test()` in the CampsPP environment, or
# `julia test/runtests.jl` from the repository root).
#
# Two suites, kept as they were when DisentangleCAMPS was a separate package: each runs in
# its own module, so their top-level definitions don't meet. DisentangleCAMPS runs first:
# its dense-reference self-test needs a module where QuantumClifford is not loaded yet.
# One wrapping testset, so a failure in one suite doesn't hide the other.
using Test

@testset "CampsPP package" begin
  @eval module DisentangleCAMPSSuite
    include($(joinpath(@__DIR__, "DisentangleCAMPS", "runtests.jl")))
  end
  @eval module CampsPPSuite
    include($(joinpath(@__DIR__, "campspp_suite.jl")))
  end
end
