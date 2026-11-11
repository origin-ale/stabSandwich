# The root project is not a package, so `Pkg.test()` cannot run here.
# This driver runs the test suites of the local packages:
#   julia test/runtests.jl
using Pkg

const ROOT = dirname(@__DIR__)
const PACKAGES = ["DisentangleCAMPS", "CampsPP"]

failed = String[]
for pkg in PACKAGES
  Pkg.activate(joinpath(ROOT, pkg))
  try
    Pkg.test()
  catch e
    e isa Pkg.Types.PkgError || rethrow()
    push!(failed, pkg)
  end
end

isempty(failed) || error("Test failures in: " * join(failed, ", "))
