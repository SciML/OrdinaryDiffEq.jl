ode_defaults = Dict(
    :dt_NaN => Verbosity.Warn(),
    :init_NaN => Verbosity.Warn(),
    :rosenbrock_no_differential_states => Verbosity.Warn(),
    :dense_output_saveat => Verbosity.Warn(),
    :alg_switch => Verbosity.Warn(),
    :mismatched_input_output_type => Verbosity.Warn(),
    :shampine_dt => Verbosity.Warn(),
    :unlimited_dt => Verbosity.Warn()
)

mutable struct ODEErrorControlVerbosity
    dt_NaN::Verbosity.Type
    init_NaN::Verbosity.Type
    dense_output_saveat::Verbosity.Type

    function ODEErrorControlVerbosity(;
            dt_NaN = ode_defaults[:dt_NaN], init_NaN = ode_defaults[:init_NaN], dense_output_saveat = ode_defaults[:dense_output_saveat])
            @info "here"
        new(dt_NaN, init_NaN, dense_output_saveat)
    end
end

function ODEErrorControlVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.Default() => ODEErrorControlVerbosity()

        Verbosity.None() => ODEErrorControlVerbosity(;NamedTuple{fieldnames(ODEErrorControlVerbosity)}(fill(
            Verbosity.None(),
            length(fieldnames(ODEErrorControlVerbosity))))...)

        Verbosity.Info() => ODEErrorControlVerbosity(;NamedTuple{fieldnames(ODEErrorControlVerbosity)}(fill(
            Verbosity.Info(),
            length(fieldnames(ODEErrorControlVerbosity))))...)

        Verbosity.Warn() => ODEErrorControlVerbosity(;NamedTuple{fieldnames(ODEErrorControlVerbosity)}(fill(
            Verbosity.Warn(),
            length(fieldnames(ODEErrorControlVerbosity))))...)

        Verbosity.Error() => ODEErrorControlVerbosity(;NamedTuple{fieldnames(ODEErrorControlVerbosity)}(fill(
            Verbosity.Error(),
            length(fieldnames(ODEErrorControlVerbosity))))...)

        Verbosity.Edge() => ODEErrorControlVerbosity()

        _ => @error "$verbose is not a valid choice for verbosity."
    end
end

mutable struct ODEPerformanceVerbosity
    alg_switch::Verbosity.Type
    mismatched_input_output_type::Verbosity.Type

    function ODEPerformanceVerbosity(;alg_switch = ode_defaults[:alg_switch],
            mismatched_input_output_type = ode_defaults[:mismatched_input_output_type])
        new(alg_switch, mismatched_input_output_type)
    end
end

function ODEPerformanceVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODEPerformanceVerbosity(;
            NamedTuple{fieldnames(ODEPerformanceVerbosity)}(fill(
                Verbosity.None(),
                length(fieldnames(ODEPerformanceVerbosity))))...)

        Verbosity.Info() => ODEPerformanceVerbosity(;
            NamedTuple{fieldnames(ODEPerformanceVerbosity)}(fill(
                Verbosity.Info(),
                length(fieldnames(ODEPerformanceVerbosity))))...)

        Verbosity.Warn() => ODEPerformanceVerbosity(;
            NamedTuple{fieldnames(ODEPerformanceVerbosity)}(fill(
                Verbosity.Warn(),
                length(fieldnames(ODEPerformanceVerbosity))))...)

        Verbosity.Error() => ODEPerformanceVerbosity(;
            NamedTuple{fieldnames(ODEPerformanceVerbosity)}(fill(
                Verbosity.Error(),
                length(fieldnames(ODEPerformanceVerbosity))))...)

        Verbosity.Default() => ODEPerformanceVerbosity()

        _ => @error "Not a valid choice for verbosity."
    end
end

mutable struct ODENumericalVerbosity
    rosenbrock_no_differential_states::Verbosity.Type
    shampine_dt::Verbosity.Type
    unlimited_dt::Verbosity.Type
    function ODENumericalVerbosity(;
            rosenbrock_no_differential_states = ode_defaults[:rosenbrock_no_differential_states],
            shampine_dt = ode_defaults[:shampine_dt],
            unlimited_dt = ode_defaults[:unlimited_dt])
        new(rosenbrock_no_differential_states, shampine_dt, unlimited_dt)
    end
end

function ODENumericalVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.None() => ODENumericalVerbosity(;
            NamedTuple{fieldnames(ODENumericalVerbosity)}(fill(
                Verbosity.None(),
                length(fieldnames(ODENumericalVerbosity))))...)

        Verbosity.Info() => OODENumericalVerbosity(;
            NamedTuple{fieldnames(ODENumericalVerbosity)}(fill(
                Verbosity.Info(),
                length(fieldnames(ODENumericalVerbosity))))...)

        Verbosity.Warn() => ODENumericalVerbosity(;
            NamedTuple{fieldnames(ODENumericalVerbosity)}(fill(
                Verbosity.Warn(),
                length(fieldnames(ODENumericalVerbosity))))...)

        Verbosity.Error() => ODENumericalVerbosity(;
            NamedTuple{fieldnames(ODENumericalVerbosity)}(fill(
                Verbosity.Error(),
                length(fieldnames(ODENumericalVerbosity))))...)

        Verbosity.Default() => ODENumericalVerbosity()

        _ => @error "Not a valid choice for verbosity."
    end
end

struct ODEVerbosity{T}
    linear_verbosity::Any
    nonlinear_verbosity::Any

    error_control::ODEErrorControlVerbosity
    performance::ODEPerformanceVerbosity
    numerical::ODENumericalVerbosity
end

function ODEVerbosity(verbose::Verbosity.Type)
    @match verbose begin
        Verbosity.Default() => ODEVerbosity{true}(
            Verbosity.Default(),
            Verbosity.Default(),
            ODEErrorControlVerbosity(Verbosity.Default()),
            ODEPerformanceVerbosity(Verbosity.Default()),
            ODENumericalVerbosity(Verbosity.Default())
        )

        Verbosity.None() => ODEVerbosity{false}(
            Verbosity.None(),
            Verbosity.None(),
            ODEErrorControlVerbosity(Verbosity.None()),
            ODEPerformanceVerbosity(Verbosity.None()),
            ODENumericalVerbosity(Verbosity.None())
        )

        Verbosity.All() => ODEVerbosity{true}(
            Verbosity.Default(),
            Verbosity.Default(),
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

    ODEVerbosity{true}(linear_verbosity, nonlinear_verbosity, error_control_verbosity,
        performance_verbosity, numerical_verbosity)
end