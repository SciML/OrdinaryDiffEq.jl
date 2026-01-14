@verbosity_specifier ODEVerbosity begin
    toggles = (
        :linear_verbosity, :nonlinear_verbosity,
        :dt_NaN, :init_NaN, :dense_output_saveat, :max_iters, :dt_min_unstable, :instability,
        :newton_convergence, :step_rejected, :step_accepted, :convergence_limit,
        :alg_switch, :stiff_detection, :mismatched_input_output_type, :jacobian_update,
        :w_factorization, :newton_iterations,
        :rosenbrock_no_differential_states, :shampine_dt, :unlimited_dt, :dt_epsilon,
        :stability_check, :near_singular,
        :sensitivity_vjp_choice,
    )

    presets = (
        None = (
            linear_verbosity = None(),
            nonlinear_verbosity = None(),
            dt_NaN = Silent(),
            init_NaN = Silent(),
            dense_output_saveat = Silent(),
            max_iters = Silent(),
            dt_min_unstable = Silent(),
            instability = Silent(),
            newton_convergence = Silent(),
            step_rejected = Silent(),
            step_accepted = Silent(),
            convergence_limit = Silent(),
            alg_switch = Silent(),
            stiff_detection = Silent(),
            mismatched_input_output_type = Silent(),
            jacobian_update = Silent(),
            w_factorization = Silent(),
            newton_iterations = Silent(),
            rosenbrock_no_differential_states = Silent(),
            shampine_dt = Silent(),
            unlimited_dt = Silent(),
            dt_epsilon = Silent(),
            stability_check = Silent(),
            near_singular = Silent(),
            sensitivity_vjp_choice = Silent(),
        ),
        Minimal = (
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
            near_singular = WarnLevel(),
            sensitivity_vjp_choice = Silent(),
        ),
        Standard = (
            linear_verbosity = Minimal(),
            nonlinear_verbosity = Minimal(),
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
            near_singular = Silent(),
            sensitivity_vjp_choice = Silent(),
        ),
        Detailed = (
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
            near_singular = WarnLevel(),
            sensitivity_vjp_choice = WarnLevel(),
        ),
        All = (
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
            near_singular = WarnLevel(),
            sensitivity_vjp_choice = WarnLevel(),
        ),
    )

    groups = (
        error_control = (
            :dt_NaN, :init_NaN, :dense_output_saveat, :max_iters, :dt_min_unstable,
            :instability, :newton_convergence, :step_rejected, :step_accepted, :convergence_limit,
        ),
        performance = (
            :alg_switch, :stiff_detection, :mismatched_input_output_type, :jacobian_update,
            :w_factorization, :newton_iterations,
        ),
        numerical = (
            :rosenbrock_no_differential_states, :shampine_dt, :unlimited_dt, :dt_epsilon,
            :stability_check, :near_singular,
        ),
        sensitivity = (
            :sensitivity_vjp_choice,
        ),
    )
end

"""
    ODEVerbosity <: AbstractVerbositySpecifier

Verbosity configuration for OrdinaryDiffEq.jl solvers, providing fine-grained control over
diagnostic messages, warnings, and errors during ODE solution.

# Fields

## Solver Verbosity
- `linear_verbosity`: Verbosity configuration for linear solvers
- `nonlinear_verbosity`: Verbosity configuration for nonlinear solvers

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

## Sensitivity Group
- `sensitivity_vjp_choice`: Messages about VJP choice in sensitivity analysis (used by SciMLSensitivity.jl)

# Constructors

    ODEVerbosity(preset::AbstractVerbosityPreset)

Create an `ODEVerbosity` using a preset configuration:
- `SciMLLogging.None()`: All messages disabled
- `SciMLLogging.Minimal()`: Only critical errors and fatal issues
- `SciMLLogging.Standard()`: Balanced verbosity (default)
- `SciMLLogging.Detailed()`: Comprehensive debugging information
- `SciMLLogging.All()`: Maximum verbosity

    ODEVerbosity(; preset=nothing, error_control=nothing, performance=nothing, numerical=nothing, kwargs...)

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
function ODEVerbosity end

const DEFAULT_VERBOSE = ODEVerbosity()

@inline function _process_verbose_param(verbose::SciMLLogging.AbstractVerbosityPreset)
    return ODEVerbosity(verbose)
end

@inline function _process_verbose_param(verbose::Bool)
    return verbose ? DEFAULT_VERBOSE : ODEVerbosity(SciMLLogging.None())
end

@inline _process_verbose_param(verbose::ODEVerbosity) = verbose
