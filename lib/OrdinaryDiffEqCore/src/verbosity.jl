ode_defaults = Dict(
    :dt_NaN => Verbosity.Warn(),
    :init_NaN => Verbosity.Warn(),
    :rosenbrock_no_differential_states => Verbosity.Warn(),
    :dense_output_saveat => Verbosity.Warn(),
    :alg_switch => Verbosity.Warn(),
    :mismatched_input_output_type => Verbosity.Warn()
)

mutable struct ODEErrorControlVerbosity
    dt_NaN::Verbosity.Type
    init_NaN::Verbosity.Type
    rosenbrock_no_differential_states::Verbosity.Type
    dense_output_saveat::Verbosity.Type

    function ODEErrorControlVerbosity(;
            dt_NaN = defaults[:dt_NaN], init_NaN = defaults[:init_NaN],
            rosenbrock_no_differential_states = defaults[:rosenbrock_no_differential_states], dense_output_saveat = defaults[:dense_output_saveat])
        new(dt_NaN, init_NaN, rosenbrock_no_differential_states, dense_output_saveat)
    end
end

function ODEErrorControlVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODEErrorControlVerbosity(fill(
            Verbosity.None(), length(fieldnames(ODEErrorControlVerbosity)))...)

        Verbosity.Info() => ODEErrorControlVerbosity(fill(
            Verbosity.Info(), length(fieldnames(ODEErrorControlVerbosity)))...)

        Verbosity.Warn() => ODEErrorControlVerbosity(fill(
            Verbosity.Warn(), length(fieldnames(ODEErrorControlVerbosity)))...)

        Verbosity.Error() => ODEErrorControlVerbosity(fill(
            Verbosity.Error(), length(fieldnames(ODEErrorControlVerbosity)))...)

        Verbosity.Default() => ODEErrorControlVerbosity()

        Verbosity.Edge() => ODEErrorControlVerbosity()

        _ => @error "Not a valid choice for verbosity."
    end
end

mutable struct ODEPerformanceVerbosity
    alg_switch::Verbosity.Type
    mismatched_input_output_type::Verbosity.Type

    function ODEPerformanceVerbosity(; alg_switch = defaults[:alg_switch],
            mismatched_input_output_type = defaults[:mismatched_input_output_type])
        new(alg_switch, mismatched_input_output_type)
    end
end

function ODEPerformanceVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODEPerformanceVerbosity(fill(
            Verbosity.None(), length(fieldnames(ODEPerformanceVerbosity)))...)

        Verbosity.Info() => ODEPerformanceVerbosity(fill(
            Verbosity.Info(), length(fieldnames(ODEPerformanceVerbosity)))...)

        Verbosity.Warn() => ODEPerformanceVerbosity(fill(
            Verbosity.Warn(), length(fieldnames(ODEPerformanceVerbosity)))...)

        Verbosity.Error() => ODEPerformanceVerbosity(fill(
            Verbosity.Error(), length(fieldnames(ODEPerformanceVerbosity)))...)

        Verbosity.Default() => ODEPerformanceVerbosity()

        _ => @error "Not a valid choice for verbosity."
    end
end

mutable struct ODENumericalVerbosity
    @add_kwonly function ODENumericalVerbosity()
        new()
    end
end

function ODENumericalVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODENumericalVerbosity(fill(
            Verbosity.None(), length(fieldnames(ODENumericalVerbosity)))...)

        Verbosity.Info() => ODENumericalVerbosity(fill(
            Verbosity.None(), length(fieldnames(ODENumericalVerbosity)))...)

        Verbosity.Warn() => ODENumericalVerbosity(fill(
            Verbosity.Warn(), length(fieldnames(ODENumericalVerbosity)))...)

        Verbosity.Error() => ODENumericalVerbosity(fill(
            Verbosity.Error(), length(fieldnames(ODENumericalVerbosity)))...)

        Verbosity.Default() => ODENumericalVerbosity()

        _ => @error "Not a valid choice for verbosity."
    end
end

struct ODEVerbosity{T} <: AbstractVerbositySpecifier{T}
    linear_verbosity
    nonlinear_verbosity

    error_control::ODEErrorControlVerbosity
    performance::ODEPerformanceVerbosity
    numerical::ODENumericalVerbosity
end

function ODEVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.Default() => ODEVerbosity{true}(
            LinearVerbosity(Verbosity.Default()),
            NonlinearVerbosity(Verbosity.Default()),
            ODEErrorControlVerbosity(Verbosity.Default()),
            ODEPerformanceVerbosity(Verbosity.Default()),
            ODENumericalVerbosity(Verbosity.Default())
        )

        Verbosity.None() => ODEVerbosity{false}(
            LinearVerbosity(Verbosity.None()),
            NonlinearVerbosity(Verbosity.None()),
            ODEErrorControlVerbosity(Verbosity.None()),
            ODEPerformanceVerbosity(Verbosity.None()),
            ODENumericalVerbosity(Verbosity.None())
        )

        Verbosity.All() => ODEVerbosity{true}(
            LinearVerbosity(Verbosity.All()),
            NonlinearVerbosity(Verbosity.All()),
            ODEErrorControlVerbosity(Verbosity.Info()),
            ODEPerformanceVerbosity(Verbosity.Info()),
            ODENumericalVerbosity(Verbosity.Info())
        )

        _ => @error "Not a valid choice for verbosity."
    end
end

function ODEVerbosity(;
        error_control = Verbosity.Default(), performance = Verbosity.Default(),
        numerical = Verbosity.Default(), linear_verbosity = Verbosity.Default(),
        nonlinear_verbosity = Verbosity.Default(), kwargs...)
    if error_control isa Verbosity.Type
        error_control_verbosity = ODEErrorControlVerbosity(error_control)
    else
        error_control_verbosity = error_control
    end

    if performance isa Verbosity.Type
        performance_verbosity = ODEPerformanceVerbosity(performance)
    else
        performance_verbosity = performance
    end

    if numerical isa Verbosity.Type
        numerical_verbosity = ODENumericalVerbosity(numerical)
    else
        numerical_verbosity = numerical
    end

    if !isempty(kwargs)
        for (key, value) in pairs(kwargs)
            if hasfield(ODEErrorControlVerbosity, key)
                setproperty!(error_control_verbosity, key, value)
            elseif hasfield(ODEPerformanceVerbosity, key)
                setproperty!(performance_verbosity, key, value)
            elseif hasfield(ODENumericalVerbosity, key)
                setproperty!(numerical_verbosity, key, value)
            else
                error("$key is not a recognized verbosity toggle.")
            end
        end
    end

    ODEVerbosity{true}(linear, nonlinear, error_control_verbosity,
        performance_verbosity, numerical_verbosity)
end