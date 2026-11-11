using DisentangleCAMPS

import Logging

struct NoDeprecationLogger <: Logging.AbstractLogger
    logger::Logging.AbstractLogger
end

Logging.min_enabled_level(l::NoDeprecationLogger) = Logging.Debug
Logging.shouldlog(l::NoDeprecationLogger, level, mod, group, id) =
    !(level == Logging.Warn && group == :depwarn)
Logging.handle_message(l::NoDeprecationLogger, args...; kwargs...) =
    Logging.handle_message(l.logger, args...; kwargs...)

Logging.global_logger(NoDeprecationLogger(Logging.global_logger()))

using Test

# One wrapping testset: a failing file no longer aborts the run and hides the rest.
@testset "DisentangleCAMPS" begin
  @testset "Dense reference (self-test)" begin
    include("TestDenseReference.jl")
  end
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
end