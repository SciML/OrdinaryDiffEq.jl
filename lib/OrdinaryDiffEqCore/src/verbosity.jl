import SciMLBase
import SciMLBase.@data
import SciMLBase.@match
using LoggingExtras

@data Verbosity begin
    None
    Edge
    Info
    Warn
    Error
    All
    Default
    Level(Int)
end


mutable struct ODEErrorControlVerbosity
    thing1::Any
    thing2::Any
end

function ODEErrorControlVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODEErrorControlVerbosity(Verbosity.None(), Verbosity.None())

        Verbosity.Warn() => ODEErrorControlVerbosity(Verbosity.Warn(), Verbosity.Warn())

        Verbosity.Error() => ODEErrorControlVerbosity(Verbosity.Error(), Verbosity.Error())

        Verbosity.Default() => ODEErrorControlVerbosity(Verbosity.Info(), Verbosity.Error())

        Verbosity.Edge() => ODEErrorControlVerbosity(Verbosity.Info(), Verbosity.Warn())

        _ => @error "Not a valid choice for verbosity."
    end
end

mutable struct ODEPerformanceVerbosity
    thing3::Any
    thing4::Any
end

function ODEPerformanceVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODEPerformanceVerbosity(Verbosity.None(), Verbosity.None())

        Verbosity.Warn() => ODEPerformanceVerbosity(Verbosity.Warn(), Verbosity.Warn())

        Verbosity.Error() => ODEPerformanceVerbosity(Verbosity.Error(), Verbosity.Error())

        Verbosity.Default() => ODEPerformanceVerbosity(Verbosity.Warn(), Verbosity.Error())

        _ => @error "Not a valid choice for verbosity."
    end
end

mutable struct ODENumericalVerbosity
    thing5::Any
    thing6::Any
end

function ODENumericalVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODENumericalVerbosity(Verbosity.None(), Verbosity.None())

        Verbosity.Warn() => ODENumericalVerbosity(Verbosity.Warn(), Verbosity.Warn())

        Verbosity.Error() => ODENumericalVerbosity(Verbosity.Error(), Verbosity.Error())

        Verbosity.Default() => ODENumericalVerbosity(Verbosity.Warn(), Verbosity.Error())

        _ => @error "Not a valid choice for verbosity."
    end
end

struct ODEVerbosity{T}
    error_control::Union{ODEErrorControlVerbosity, Nothing}
    performance::Union{ODEPerformanceVerbosity, Nothing}
    numerical::Union{ODENumericalVerbosity, Nothing}
end

function ODEVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.Default() => ODEVerbosity{true}(
            ODEErrorControlVerbosity(Verbosity.Default()),
            ODEPerformanceVerbosity(Verbosity.Default()),
            ODENumericalVerbosity(Verbosity.Default())
        )

        Verbosity.None() => ODEVerbosity{false}(nothing, nothing, nothing)

        Verbosity.All() => ODEVerbosity{true}(
            ODEErrorControlVerbosity(Verbosity.All()),
            ODEPerformanceVerbosity(Verbosity.All()),
            ODENumericalVerbosity(Verbosity.All())
        )

        _ => @error "Not a valid choice for verbosity."
    end
end

function (verbose::ODEVerbosity{true})(message, option, group, values...)
    level = get_message_level(verbose, option, group)

    if !isnothing(level)
        Base.@logmsg level message values
    end
end

function (verbose::ODEVerbosity{true})(f::Function, option, group)
    level = get_message_level(verbose, option, group)

    if !isnothing(level)
        message = f()
        Base.@logmsg level message
    end
end

function (verbose::ODEVerbosity{false})(f::Function, option, group)
end

function (verbose::ODEVerbosity{false})(message, option, group)
end

function get_message_level(verbose::ODEVerbosity{true}, option, group)
    group = getproperty(verbose, group)
    opt_level = getproperty(group, option)

    @match opt_level begin
        Verbosity.None() => nothing
        Verbosity.Info() => Logging.Info
        Verbosity.Warn() => Logging.Warn
        Verbosity.Error() => Logging.Error
        Verbosity.Level(i) => Logging.LogLevel(i)
    end
end

function ODELogger(; info_repl = true, warn_repl = true, error_repl = true,
        info_file = nothing, warn_file = nothing, error_file = nothing)
    info_sink = isnothing(info_file) ? NullLogger() : FileLogger(info_file)
    warn_sink = isnothing(warn_file) ? NullLogger() : FileLogger(warn_file)
    error_sink = isnothing(error_file) ? NullLogger() : FileLogger(error_file)

    repl_filter = EarlyFilteredLogger(current_logger()) do log
        if log.level == Logging.Info && info_repl
            return true
        end

        if log.level == Logging.Warn && warn_repl
            return true
        end

        if log.level == Logging.Error && error_repl
            return true
        end

        return false
    end

    info_filter = EarlyFilteredLogger(info_sink) do log
        log.level == Logging.Info
    end

    warn_filter = EarlyFilteredLogger(warn_sink) do log
        log.level == Logging.Warn
    end

    error_filter = EarlyFilteredLogger(error_sink) do log
        log.level == Logging.Error
    end

    TeeLogger(repl_filter, info_filter, warn_filter, error_filter)
end


