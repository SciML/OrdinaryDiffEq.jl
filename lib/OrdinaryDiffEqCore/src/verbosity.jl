"""
    ODEVerbosity <: AbstractVerbositySpecifier

Verbosity configuration for OrdinaryDiffEq.jl solvers, providing fine-grained control over
diagnostic messages, warnings, and errors during ODE solution.

# Fields

## Error Control Group
- `dt_NaN`: Messages when time step becomes NaN
- `init_NaN`: Messages when initial conditions contain NaN
- `dense_output_saveat`: Messages about dense output with saveat
- `max_iters`: Messages when maximum iterations are reached
- `dt_min_unstable`: Messages when time step becomes too small/unstable
- `instability`: Messages when numerical instability is detected
- `newton_convergence`: Messages when Newton iteration fails to converge
- `step_rejected`: Messages when adaptive steps are rejected
- `step_accepted`: Messages when adaptive steps are accepted
- `convergence_limit`: Messages when convergence at floating point precision limit

## Performance Group
- `alg_switch`: Messages when algorithm switching occurs
- `stiff_detection`: Messages when stiffness is detected
- `mismatched_input_output_type`: Messages when input/output types don't match
- `jacobian_update`: Messages when Jacobian matrix is computed/updated
- `w_factorization`: Messages when W matrix is factorized
- `newton_iterations`: Messages about Newton iteration progress

## Numerical Group
- `rosenbrock_no_differential_states`: Messages when Rosenbrock has no differential states
- `shampine_dt`: Messages about Shampine time step selection
- `unlimited_dt`: Messages when time step is unlimited
- `dt_epsilon`: Messages when timestep goes below floating point epsilon
- `stability_check`: Messages about stability checks in extrapolation methods
- `near_singular`: Messages when Jacobian/mass matrix appears near-singular

## Solver Verbosity Groups
- `linear_verbosity`: Verbosity configuration for linear solvers
- `nonlinear_verbosity`: Verbosity configuration for nonlinear solvers

# Constructors

    ODEVerbosity(preset::AbstractVerbosityPreset)

Create an `ODEVerbosity` using a preset configuration:
- `SciMLLogging.None()`: All messages disabled
- `SciMLLogging.Minimal()`: Only critical errors and fatal issues
- `SciMLLogging.Standard()`: Balanced verbosity (default)
- `SciMLLogging.Detailed()`: Comprehensive debugging information
- `SciMLLogging.All()`: Maximum verbosity

    ODEVerbosity(; error_control=nothing, performance=nothing, numerical=nothing, linear_verbosity=nothing, nonlinear_verbosity=nothing, kwargs...)

Create an `ODEVerbosity` with group-level or individual field control.

# Examples

```julia
# Use a preset
verbose = ODEVerbosity(SciMLLogging.Standard())

# Set entire groups
verbose = ODEVerbosity(
    error_control = SciMLLogging.WarnLevel(),
    numerical = SciMLLogging.InfoLevel()
)

# Set individual fields
verbose = ODEVerbosity(
    dt_NaN = SciMLLogging.ErrorLevel(),
    alg_switch = SciMLLogging.InfoLevel()
)

# Mix group and individual settings
verbose = ODEVerbosity(
    numerical = SciMLLogging.InfoLevel(),  # Set all numerical to InfoLevel
    unlimited_dt = SciMLLogging.ErrorLevel()  # Override specific field
)
```
"""
@concrete struct ODEVerbosity <: AbstractVerbositySpecifier
    # Solver verbosity
    linear_verbosity
    nonlinear_verbosity
    # Error control
    dt_NaN
    init_NaN
    dense_output_saveat
    max_iters
    dt_min_unstable
    instability
    newton_convergence
    step_rejected
    step_accepted
    convergence_limit
    # Performance
    alg_switch
    stiff_detection
    mismatched_input_output_type
    jacobian_update
    w_factorization
    newton_iterations
    # Numerical
    rosenbrock_no_differential_states
    shampine_dt
    unlimited_dt
    dt_epsilon
    stability_check
    near_singular
end

# Group classifications
const error_control_options = (:dt_NaN, :init_NaN, :dense_output_saveat, :max_iters, :dt_min_unstable, :instability, :newton_convergence, :step_rejected, :step_accepted, :convergence_limit)
const performance_options = (:alg_switch, :stiff_detection, :mismatched_input_output_type, :jacobian_update, :w_factorization, :newton_iterations)
const numerical_options = (:rosenbrock_no_differential_states, :shampine_dt, :unlimited_dt, :dt_epsilon, :stability_check, :near_singular)

# Runtime helper for complex verbosity construction path
function _build_ode_verbosity_runtime(
    error_control, performance, numerical,
    linear_verbosity, nonlinear_verbosity, kwargs
)
    # Validate group arguments
    if error_control !== nothing && !(error_control isa AbstractMessageLevel)
        throw(ArgumentError("error_control must be a SciMLLogging.AbstractMessageLevel, got $(typeof(error_control))"))
    end
    if performance !== nothing && !(performance isa AbstractMessageLevel)
        throw(ArgumentError("performance must be a SciMLLogging.AbstractMessageLevel, got $(typeof(performance))"))
    end
    if numerical !== nothing && !(numerical isa AbstractMessageLevel)
        throw(ArgumentError("numerical must be a SciMLLogging.AbstractMessageLevel, got $(typeof(numerical))"))
    end

    # Validate individual kwargs
    for (key, value) in pairs(kwargs)
        if !(key in error_control_options || key in performance_options ||
             key in numerical_options)
            throw(ArgumentError("Unknown verbosity option: $key. Valid options are: $(tuple(error_control_options..., performance_options..., numerical_options...))"))
        end
        if !(value isa AbstractMessageLevel)
            throw(ArgumentError("$key must be a SciMLLogging.AbstractMessageLevel, got $(typeof(value))"))
        end
    end

    # Build arguments using NamedTuple for type stability
    default_args = (
        linear_verbosity = linear_verbosity === nothing ? Minimal() : linear_verbosity,
        nonlinear_verbosity = nonlinear_verbosity === nothing ? Minimal() : nonlinear_verbosity,
        dt_NaN = WarnLevel(),
        init_NaN = WarnLevel(),
        dense_output_saveat = WarnLevel(),
        max_iters = WarnLevel(),
        dt_min_unstable = WarnLevel(),
        instability = WarnLevel(),
        newton_convergence = Silent(),
        step_rejected = Silent(),
        step_accepted = Silent(),
        convergence_limit = Silent(),
        alg_switch = Silent(),
        stiff_detection = Silent(),
        mismatched_input_output_type = WarnLevel(),
        jacobian_update = Silent(),
        w_factorization = Silent(),
        newton_iterations = Silent(),
        rosenbrock_no_differential_states = WarnLevel(),
        shampine_dt = Silent(),
        unlimited_dt = WarnLevel(),
        dt_epsilon = Silent(),
        stability_check = Silent(),
        near_singular = Silent()
    )

    # Apply group-level settings - done explicitly for type stability
    final_args = if error_control !== nothing || performance !== nothing ||
                    numerical !== nothing
        (
            linear_verbosity = default_args.linear_verbosity,
            nonlinear_verbosity = default_args.nonlinear_verbosity,
            # Error control group
            dt_NaN = error_control !== nothing ? error_control : default_args.dt_NaN,
            init_NaN = error_control !== nothing ? error_control : default_args.init_NaN,
            dense_output_saveat = error_control !== nothing ? error_control : default_args.dense_output_saveat,
            max_iters = error_control !== nothing ? error_control : default_args.max_iters,
            dt_min_unstable = error_control !== nothing ? error_control : default_args.dt_min_unstable,
            instability = error_control !== nothing ? error_control : default_args.instability,
            newton_convergence = error_control !== nothing ? error_control : default_args.newton_convergence,
            step_rejected = error_control !== nothing ? error_control : default_args.step_rejected,
            step_accepted = error_control !== nothing ? error_control : default_args.step_accepted,
            convergence_limit = error_control !== nothing ? error_control : default_args.convergence_limit,
            # Performance group
            alg_switch = performance !== nothing ? performance : default_args.alg_switch,
            stiff_detection = performance !== nothing ? performance : default_args.stiff_detection,
            mismatched_input_output_type = performance !== nothing ? performance : default_args.mismatched_input_output_type,
            jacobian_update = performance !== nothing ? performance : default_args.jacobian_update,
            w_factorization = performance !== nothing ? performance : default_args.w_factorization,
            newton_iterations = performance !== nothing ? performance : default_args.newton_iterations,
            # Numerical group
            rosenbrock_no_differential_states = numerical !== nothing ? numerical : default_args.rosenbrock_no_differential_states,
            shampine_dt = numerical !== nothing ? numerical : default_args.shampine_dt,
            unlimited_dt = numerical !== nothing ? numerical : default_args.unlimited_dt,
            dt_epsilon = numerical !== nothing ? numerical : default_args.dt_epsilon,
            stability_check = numerical !== nothing ? numerical : default_args.stability_check,
            near_singular = numerical !== nothing ? numerical : default_args.near_singular
        )
    else
        default_args
    end

    # Apply individual overrides
    if !isempty(kwargs)
        final_args = merge(final_args, NamedTuple(kwargs))
    end

    ODEVerbosity(values(final_args)...)
end

# Optionally-generated function for ODEVerbosity construction
# Uses compile-time type information to optimize the common default case
function _build_ode_verbosity(
    error_control,
    performance,
    numerical,
    linear_verbosity,
    nonlinear_verbosity,
    kwargs
)
    if @generated
        # Generated path: In this block, we're at compile time
        # error_control, performance, numerical, etc. are the actual type values
        # Check if all group params are Nothing and kwargs is empty (fast default path)
        if error_control === Nothing && performance === Nothing && numerical === Nothing &&
           linear_verbosity === Nothing && nonlinear_verbosity === Nothing &&
           kwargs <: NamedTuple{()}
            # Return an expression that constructs the default directly
            return quote
                ODEVerbosity(
                    Minimal(),          # linear_verbosity
                    Minimal(),          # nonlinear_verbosity
                    WarnLevel(),        # dt_NaN
                    WarnLevel(),        # init_NaN
                    WarnLevel(),        # dense_output_saveat
                    WarnLevel(),        # max_iters
                    WarnLevel(),        # dt_min_unstable
                    WarnLevel(),        # instability
                    Silent(),           # newton_convergence
                    Silent(),           # step_rejected
                    Silent(),           # step_accepted
                    Silent(),           # convergence_limit
                    Silent(),           # alg_switch
                    Silent(),           # stiff_detection
                    WarnLevel(),        # mismatched_input_output_type
                    Silent(),           # jacobian_update
                    Silent(),           # w_factorization
                    Silent(),           # newton_iterations
                    WarnLevel(),        # rosenbrock_no_differential_states
                    Silent(),           # shampine_dt
                    WarnLevel(),        # unlimited_dt
                    Silent(),           # dt_epsilon
                    Silent(),           # stability_check
                    Silent()            # near_singular
                )
            end
        else
            # For non-default cases, delegate to runtime logic
            return quote
                _build_ode_verbosity_runtime(
                    error_control, performance, numerical,
                    linear_verbosity, nonlinear_verbosity, kwargs
                )
            end
        end
    else
        # Runtime fallback: error_control, performance, etc. are VALUES
        if error_control === nothing && performance === nothing && numerical === nothing &&
           linear_verbosity === nothing && nonlinear_verbosity === nothing &&
           isempty(kwargs)
            # Fast default path at runtime
            ODEVerbosity(
                Minimal(),          # linear_verbosity
                Minimal(),          # nonlinear_verbosity
                WarnLevel(),        # dt_NaN
                WarnLevel(),        # init_NaN
                WarnLevel(),        # dense_output_saveat
                WarnLevel(),        # max_iters
                WarnLevel(),        # dt_min_unstable
                WarnLevel(),        # instability
                Silent(),           # newton_convergence
                Silent(),           # step_rejected
                Silent(),           # step_accepted
                Silent(),           # convergence_limit
                Silent(),           # alg_switch
                Silent(),           # stiff_detection
                WarnLevel(),        # mismatched_input_output_type
                Silent(),           # jacobian_update
                Silent(),           # w_factorization
                Silent(),           # newton_iterations
                WarnLevel(),        # rosenbrock_no_differential_states
                Silent(),           # shampine_dt
                WarnLevel(),        # unlimited_dt
                Silent(),           # dt_epsilon
                Silent(),           # stability_check
                Silent()            # near_singular
            )
        else
            # Complex path
            _build_ode_verbosity_runtime(
                error_control, performance, numerical,
                linear_verbosity, nonlinear_verbosity, kwargs
            )
        end
    end
end

function ODEVerbosity(;
        error_control = nothing, performance = nothing, numerical = nothing,
        linear_verbosity = nothing, nonlinear_verbosity = nothing, kwargs...)
    _build_ode_verbosity(
        error_control, performance, numerical,
        linear_verbosity, nonlinear_verbosity,
        NamedTuple(kwargs)
    )
end

# Constructor for verbosity presets following the hierarchical levels:
# None < Minimal < Standard < Detailed < All
# Each level includes all messages from levels below it plus additional ones
function ODEVerbosity(verbose::AbstractVerbosityPreset)
    if verbose isa Minimal
        # Minimal: Only fatal errors and critical warnings
        return ODEVerbosity(
            linear_verbosity = Minimal(),
            nonlinear_verbosity = Minimal(),
            dt_NaN = WarnLevel(),
            init_NaN = WarnLevel(),
            dense_output_saveat = Silent(),
            max_iters = WarnLevel(),
            dt_min_unstable = WarnLevel(),
            instability = WarnLevel(),
            newton_convergence = WarnLevel(),
            step_rejected = Silent(),
            step_accepted = Silent(),
            convergence_limit = Silent(),
            alg_switch = Silent(),
            stiff_detection = Silent(),
            mismatched_input_output_type = Silent(),
            jacobian_update = Silent(),
            w_factorization = Silent(),
            newton_iterations = Silent(),
            rosenbrock_no_differential_states = WarnLevel(),
            shampine_dt = Silent(),
            unlimited_dt = WarnLevel(),
            dt_epsilon = Silent(),
            stability_check = Silent(),
            near_singular = WarnLevel()
        )
    elseif verbose isa Standard
        # Standard: Everything from Minimal + non-fatal warnings
        return ODEVerbosity()
    elseif verbose isa Detailed
        # Detailed: Everything from Standard + debugging/solver behavior
        return ODEVerbosity(
            linear_verbosity = Detailed(),
            nonlinear_verbosity = Detailed(),
            dt_NaN = WarnLevel(),
            init_NaN = WarnLevel(),
            dense_output_saveat = InfoLevel(),
            max_iters = WarnLevel(),
            dt_min_unstable = WarnLevel(),
            instability = WarnLevel(),
            newton_convergence = WarnLevel(),
            step_rejected = Silent(),
            step_accepted = Silent(),
            convergence_limit = InfoLevel(),
            alg_switch = InfoLevel(),
            stiff_detection = Silent(),
            mismatched_input_output_type = WarnLevel(),
            jacobian_update = InfoLevel(),
            w_factorization = InfoLevel(),
            newton_iterations = InfoLevel(),
            rosenbrock_no_differential_states = WarnLevel(),
            shampine_dt = InfoLevel(),
            unlimited_dt = WarnLevel(),
            dt_epsilon = InfoLevel(),
            stability_check = InfoLevel(),
            near_singular = WarnLevel()
        )
    elseif verbose isa All
        # All: Maximum verbosity - every possible logging message at InfoLevel
        return ODEVerbosity(
            linear_verbosity = All(),
            nonlinear_verbosity = All(),
            dt_NaN = WarnLevel(),
            init_NaN = WarnLevel(),
            dense_output_saveat = InfoLevel(),
            max_iters = WarnLevel(),
            dt_min_unstable = WarnLevel(),
            instability = WarnLevel(),
            newton_convergence = WarnLevel(),
            step_rejected = InfoLevel(),
            step_accepted = InfoLevel(),
            convergence_limit = InfoLevel(),
            alg_switch = InfoLevel(),
            stiff_detection = InfoLevel(),
            mismatched_input_output_type = InfoLevel(),
            jacobian_update = InfoLevel(),
            w_factorization = InfoLevel(),
            newton_iterations = InfoLevel(),
            rosenbrock_no_differential_states = WarnLevel(),
            shampine_dt = InfoLevel(),
            unlimited_dt = WarnLevel(),
            dt_epsilon = InfoLevel(),
            stability_check = InfoLevel(),
            near_singular = WarnLevel()
        )
    end
end

@inline function ODEVerbosity(verbose::None)
    ODEVerbosity(
        None(),
        None(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent(),
        Silent()
    )
end

const DEFAULT_VERBOSE = ODEVerbosity()

@inline function _process_verbose_param(verbose::SciMLLogging.AbstractVerbosityPreset)
    ODEVerbosity(verbose)
end

@inline function _process_verbose_param(verbose::Bool)
    verbose ? DEFAULT_VERBOSE : ODEVerbosity(SciMLLogging.None())
end 

@inline _process_verbose_param(verbose::ODEVerbosity) = verbose