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

include("TestDisentangler.jl")
include("TestCompatibility.jl")
# include("TestEvolve.jl")