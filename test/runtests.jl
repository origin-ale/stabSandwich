# The root project is not a package, so `Pkg.test()` cannot run here.
# This driver runs the test suites of the local packages, then the tests of the root-level
# script helpers (scripts/, not part of any package) in the root environment:
#   julia test/runtests.jl
using Pkg

const ROOT = dirname(@__DIR__)
const PACKAGES = ["CampsPP"]   # DisentangleCAMPS was merged into CampsPP
# Run in a child process each, in the root environment, with 2 threads (setup_threads)
const ROOT_TESTS = ["TestScriptSetup.jl"]

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

for file in ROOT_TESTS
  cmd = `$(Base.julia_cmd()) --project=$ROOT --startup-file=no -t 2 $(joinpath(@__DIR__, file))`
  success(pipeline(ignorestatus(cmd); stdout, stderr)) || push!(failed, file)
end

isempty(failed) || error("Test failures in: " * join(failed, ", "))
