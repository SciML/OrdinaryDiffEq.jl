struct EvalFunc{F} <: Function
    f::F
end
(f::EvalFunc)(args...) = f.f(args...)

NO_TSPAN_PROBS = Union{
    AbstractLinearProblem, AbstractNonlinearProblem,
    AbstractIntegralProblem, AbstractSteadyStateProblem,
    AbstractJumpProblem,
}

"""
    has_callbacks(kwargs)

Check if there are any callbacks in the kwargs. Returns `true` if callbacks are present.
"""
function has_callbacks(kwargs)
    cb = get(kwargs, :callback, nothing)
    cb === nothing && return false
    cb isa CallbackSet && return !isempty(cb)
    return true
end

"""
    merge_problem_kwargs(prob; merge_callbacks=true, kwargs…)

Merges kwargs stored in `prob.kwargs` with the provided `kwargs`, following
DiffEq's standard merging rules:
- Problem kwargs are merged first
- Passed kwargs take precedence (i.e., they override problem kwargs)
- If `merge_callbacks=true` and both prob and kwargs have callbacks, they are
  merged into a `CallbackSet` rather than one overriding the other

Returns the merged kwargs as a Base.pairs.

This function is intended for use by problem types that override `__solve` or `__init`
and need to manually handle kwargs merging that would normally be done by `solve_call`
or `init_call`.
"""
function merge_problem_kwargs(prob; merge_callbacks = true, kwargs...)

    # Special handling for callback merging
    if has_kwargs(prob)
        if merge_callbacks && haskey(prob.kwargs, :callback) && haskey(kwargs, :callback)
            kwargs_temp = NamedTuple{
                Base.diff_names(
                    Base._nt_names(values(kwargs)),
                    (:callback,)
                ),
            }(values(kwargs))
            callbacks = NamedTuple{(:callback,)}(
                (
                    DiffEqBase.CallbackSet(
                        prob.kwargs[:callback],
                        values(kwargs).callback
                    ),
                )
            )
            kwargs = merge(kwargs_temp, callbacks)
        end
        kwargs = isempty(prob.kwargs) ? kwargs : merge(values(prob.kwargs), kwargs)
    end

    return kwargs
end

function init_call(
        _prob, args...; merge_callbacks = true, kwargshandle = nothing,
        kwargs...
    )
    kwargshandle = kwargshandle === nothing ? SciMLBase.KeywordArgError : kwargshandle
    kwargshandle = has_kwargs(_prob) && haskey(_prob.kwargs, :kwargshandle) ?
        _prob.kwargs[:kwargshandle] : kwargshandle

    # Merge problem kwargs with passed kwargs
    kwargs = merge_problem_kwargs(_prob; merge_callbacks, kwargs...)

    checkkwargs(kwargshandle; kwargs...)

    return if _prob isa Union{ODEProblem, DAEProblem} && isnothing(_prob.u0) &&
            !has_callbacks(kwargs)
        build_null_integrator(_prob, args...; kwargs...)
    elseif hasfield(typeof(_prob), :f) && hasfield(typeof(_prob.f), :f) &&
            _prob.f.f isa EvalFunc
        Base.invokelatest(__init, _prob, args...; kwargs...) #::T
    else
        __init(_prob, args...; kwargs...) #::T
    end
end

function init(
        prob::AbstractDEProblem, args...; sensealg = nothing,
        u0 = nothing, p = nothing, kwargs...
    )
    if sensealg === nothing && has_kwargs(prob) && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    u0 = u0 !== nothing ? u0 : prob.u0
    p = p !== nothing ? p : prob.p

    return init_up(prob, sensealg, u0, p, args...; kwargs...)
end

function init(prob::AbstractJumpProblem, args...; kwargs...)
    return init_call(prob, args...; kwargs...)
end

function init_up(prob::AbstractDEProblem, sensealg, u0, p, args...; kwargs...)
    alg = extract_alg(args, kwargs, has_kwargs(prob) ? prob.kwargs : kwargs)
    return if isnothing(alg) || !(alg isa AbstractDEAlgorithm) # Default algorithm handling
        _prob = get_concrete_problem(
            prob, !(prob isa DiscreteProblem); alg = alg, u0 = u0,
            p = p, kwargs...
        )
        init_call(_prob, args...; kwargs...)
    else
        tstops = get(kwargs, :tstops, nothing)
        if tstops === nothing && has_kwargs(prob)
            tstops = get(prob.kwargs, :tstops, nothing)
        end
        if !(tstops isa Union{Nothing, AbstractArray, Tuple, Real}) &&
                !SciMLBase.allows_late_binding_tstops(alg)
            throw(LateBindingTstopsNotSupportedError())
        end
        _prob = get_concrete_problem(prob, isadaptive(alg); alg = alg, u0 = u0, p = p, kwargs...)
        _alg = prepare_alg(alg, _prob.u0, _prob.p, _prob)
        check_prob_alg_pairing(_prob, alg) # alg for improved inference
        if length(args) > 1
            init_call(_prob, _alg, Base.tail(args)...; kwargs...)
        else
            init_call(_prob, _alg; kwargs...)
        end
    end
end

function solve_call(
        _prob, args...; merge_callbacks = true, kwargshandle = nothing,
        kwargs...
    )
    kwargshandle = kwargshandle === nothing ? SciMLBase.KeywordArgError : kwargshandle
    kwargshandle = has_kwargs(_prob) && haskey(_prob.kwargs, :kwargshandle) ?
        _prob.kwargs[:kwargshandle] : kwargshandle

    # Merge problem kwargs with passed kwargs
    kwargs = merge_problem_kwargs(_prob; merge_callbacks, kwargs...)

    checkkwargs(kwargshandle; kwargs...)
    if isdefined(_prob, :u0)
        if _prob.u0 isa Array
            if !isconcretetype(RecursiveArrayTools.recursive_unitless_eltype(_prob.u0))
                throw(NonConcreteEltypeError(RecursiveArrayTools.recursive_unitless_eltype(_prob.u0)))
            end

            if !(eltype(_prob.u0) <: Number) && !(eltype(_prob.u0) <: Enum) &&
                    !(_prob.u0 isa AbstractVector{<:AbstractArray} && _prob isa BVProblem)
                # Allow Enums for FunctionMaps, make into a trait in the future
                # BVPs use Vector of Arrays for initial guesses
                throw(NonNumberEltypeError(eltype(_prob.u0)))
            end
        end

        if _prob.u0 === nothing && !has_callbacks(kwargs)
            return build_null_solution(_prob, args...; kwargs...)
        end
    end

    return if hasfield(typeof(_prob), :f) && hasfield(typeof(_prob.f), :f) &&
            _prob.f.f isa EvalFunc
        Base.invokelatest(__solve, _prob, args...; kwargs...) #::T
    else
        __solve(_prob, args...; kwargs...) #::T
    end
end

mutable struct NullODEIntegrator{IIP, ProbType, T, SolType, F, P} <:
    AbstractODEIntegrator{Nothing, IIP, Nothing, T}
    du::Vector{Float64}
    u::Vector{Float64}
    t::T
    prob::ProbType
    sol::SolType
    f::F
    p::P
end
function build_null_integrator(
        prob::AbstractDEProblem, args...;
        kwargs...
    )
    sol = solve(prob, args...; kwargs...)
    # The DAE initialization in `build_null_solution` may change the parameter
    # object `prob.p` via `@set!`, hence use the "new" prob instead of the "old" one.
    prob = sol.prob
    return NullODEIntegrator{
        isinplace(prob), typeof(prob), eltype(prob.tspan), typeof(sol),
        typeof(prob.f), typeof(prob.p),
    }(
        Float64[],
        Float64[],
        prob.tspan[1],
        prob,
        sol,
        prob.f,
        prob.p
    )
end
function solve!(integ::NullODEIntegrator)
    integ.t = integ.sol.t[end]
    return nothing
end
function step!(integ::NullODEIntegrator, dt = nothing, stop_at_tdt = false)
    if !isnothing(dt)
        integ.t += dt
    else
        integ.t = integ.sol.t[end]
    end
    return nothing
end
function SciMLBase.u_modified!(integ::NullODEIntegrator, u) end
SciMLBase.check_error(integ::NullODEIntegrator) = integ.sol.retcode

function hack_null_solution_init(prob)
    if SciMLBase.has_initialization_data(prob.f)
        initializeprob = prob.f.initialization_data.initializeprob
        nlsol = solve(initializeprob)
        success = SciMLBase.successful_retcode(nlsol)
        if prob.f.initialization_data.initializeprobpmap !== nothing
            @set! prob.p = prob.f.initializeprobpmap(prob, nlsol)
        end
    else
        success = true
    end
    return prob, success
end

function build_null_solution(
        prob::AbstractDEProblem, args...;
        saveat = (),
        save_everystep = true,
        save_on = true,
        save_start = save_everystep || isempty(saveat) ||
            saveat isa Number || prob.tspan[1] in saveat,
        save_end = true,
        kwargs...
    )
    ts = if saveat === ()
        if save_start && save_end
            [prob.tspan[1], prob.tspan[2]]
        elseif save_start && !save_end
            [prob.tspan[1]]
        elseif !save_start && save_end
            [prob.tspan[2]]
        else
            eltype(prob.tspan)[]
        end
    elseif saveat isa Number
        prob.tspan[1]:saveat:prob.tspan[2]
    else
        saveat
    end

    timeseries = [Float64[] for i in 1:length(ts)]

    prob, success = hack_null_solution_init(prob)
    retcode = success ? ReturnCode.Success : ReturnCode.InitialFailure
    return build_solution(prob, nothing, ts, timeseries; dense = true, retcode)
end

"""
```julia
solve(prob::AbstractDEProblem, alg::Union{AbstractDEAlgorithm,Nothing}; kwargs...)
```

## Arguments

The only positional argument is `alg` which is optional. By default, `alg = nothing`.
If `alg = nothing`, then `solve` dispatches to the DifferentialEquations.jl automated
algorithm selection (if `using DifferentialEquations` was done, otherwise it will
error with a `MethodError`).

## Keyword Arguments

The DifferentialEquations.jl universe has a large set of common arguments available
for the `solve` function. These arguments apply to `solve` on any problem type and
are only limited by limitations of the specific implementations.

Many of the defaults depend on the algorithm or the package the algorithm derives
from. Not all of the interface is provided by every algorithm.
For more detailed information on the defaults and the available options
for specific algorithms / packages, see the manual pages for the solvers of specific
problems. To see whether a specific package is compatible with the use of a
given option, see the [Solver Compatibility Chart](https://docs.sciml.ai/DiffEqDocs/stable/basics/compatibility_chart/#Solver-Compatibility-Chart)

### Default Algorithm Hinting

To help choose the default algorithm, the keyword argument `alg_hints` is
provided to `solve`. `alg_hints` is a `Vector{Symbol}` which describe the
problem at a high level to the solver. The options are:

* `:auto` vs `:nonstiff` vs `:stiff` - Denotes the equation as nonstiff/stiff.
  `:auto` allow the default handling algorithm to choose stiffness detection
  algorithms. The default handling defaults to using `:auto`.

Currently unused options include:

* `:interpolant` - Denotes that a high-precision interpolation is important.
* `:memorybound` - Denotes that the solver will be memory bound.

This functionality is derived via the benchmarks in
[SciMLBenchmarks.jl](https://github.com/SciML/SciMLBenchmarks.jl)

#### SDE Specific Alghints

* `:additive` - Denotes that the underlying SDE has additive noise.
* `:stratonovich` - Denotes that the solution should adhere to the Stratonovich
  interpretation.

### Output Control

These arguments control the output behavior of the solvers. It defaults to maximum
output to give the best interactive user experience, but can be reduced all the
way to only saving the solution at the final timepoint.

The following options are all related to output control. See the "Examples"
section at the end of this page for some example usage.

* `dense`: Denotes whether to save the extra pieces required for dense (continuous)
  output. Default is `save_everystep && isempty(saveat)` for algorithms which have
  the ability to produce dense output, i.e. by default it's `true` unless the user
  has turned off saving on steps or has chosen a `saveat` value. If `dense=false`,
  the solution still acts like a function, and `sol(t)` is a linear interpolation
  between the saved time points.
* `saveat`: Denotes specific times to save the solution at, during the solving
  phase. The solver will save at each of the timepoints in this array in the
  most efficient manner available to the solver. If only `saveat` is given, then
  the arguments `save_everystep` and `dense` are `false` by default.
  If `saveat` is given a number, then it will automatically expand to
  `tspan[1]:saveat:tspan[2]`. For methods where interpolation is not possible,
  `saveat` may be equivalent to `tstops`. The default value is `[]`.
* `save_idxs`: Denotes the indices for the components of the equation to save.
  Defaults to saving all indices. For example, if you are solving a 3-dimensional ODE,
  and given `save_idxs = [1, 3]`, only the first and third components of the
  solution will be outputted.
  Notice that of course in this case the outputted solution will be two-dimensional.
* `tstops`: Denotes *extra* times that the timestepping algorithm must step to.
  This should be used to help the solver deal with discontinuities and
  singularities, since stepping exactly at the time of the discontinuity will
  improve accuracy. If a method cannot change timesteps (fixed timestep
  multistep methods), then `tstops` will use an interpolation,
  matching the behavior of `saveat`. If a method cannot change timesteps and
  also cannot interpolate, then `tstops` must be a multiple of `dt` or else an
  error will be thrown. `tstops` may also be a function `tstops(p, tspan)`, accepting the parameter
  object and `tspan`, returning the vector of time points to stop at. Default is `[]`.
* `d_discontinuities:` Denotes locations of discontinuities in low order derivatives.
  This will force FSAL algorithms which assume derivative continuity to re-evaluate
  the derivatives at the point of discontinuity. The default is `[]`.
* `save_everystep`: Saves the result at every step.
  Default is true if `isempty(saveat)`.
* `save_on`: Denotes whether intermediate solutions are saved. This overrides the
  settings of `dense`, `saveat` and `save_everystep` and is used by some applications
  to manually turn off saving temporarily. Everyday use of the solvers should leave
  this unchanged. Defaults to `true`.
* `save_start`: Denotes whether the initial condition should be included in the solution
  type as the first timepoint. This setting overrides `saveat` when set to `false`.
  Defaults to
  `save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[1] in saveat`.
* `save_end`: Denotes whether the final condition should be included in the solution type
  as the final timepoint. This setting is overridden by other saving settings when set to
  `false`. Defaults to
  `save_everystep || isempty(saveat) || saveat isa Number || prob.tspan[2] in saveat`.
* `initialize_save`: Denotes whether to save after the callback initialization
  phase (when `u_modified=true`). Defaults to `true`.

Note that `dense` requires `save_everystep=true` and `saveat=false`. If you need
additional saving while keeping dense output, see
[the SavingCallback in the Callback Library](https://docs.sciml.ai/DiffEqCallbacks/stable/output_saving/#DiffEqCallbacks.SavingCallback).

### Stepsize Control

These arguments control the timestepping routines.

#### Basic Stepsize Control

These are the standard options for controlling stepping behavior. Error estimates
do the comparison

```math
err_{scaled} = err/(abstol + max(uprev,u)*reltol)
```

The scaled error is guaranteed to be `<1` for a given local error estimate
(note: error estimates are local unless the method specifies otherwise). `abstol`
controls the non-scaling error and thus can be thought of as the error around zero.
`reltol` scales with the size of the dependent variables and so one can interpret
`reltol=1e-3` as roughly being (locally) correct to 3 digits. Note tolerances can
be specified element-wise by passing a vector whose size matches `u0`.

* `adaptive`: Turns on adaptive timestepping for appropriate methods. Default
  is true.
* `abstol`: Absolute tolerance in adaptive timestepping. This is the tolerance
  on local error estimates, not necessarily the global error (though these quantities
  are related). Defaults to `1e-6` on deterministic equations (ODEs/DDEs/DAEs) and `1e-2`
  on stochastic equations (SDEs/RODEs).
* `reltol`: Relative tolerance in adaptive timestepping.  This is the tolerance
  on local error estimates, not necessarily the global error (though these quantities
  are related). Defaults to `1e-3` on deterministic equations (ODEs/DDEs/DAEs) and `1e-2`
  on stochastic equations (SDEs/RODEs).
* `dt`: Sets the initial stepsize. This is also the stepsize for fixed
  timestep methods. Defaults to an automatic choice if the method is adaptive.
* `dtmax`: Maximum dt for adaptive timestepping. Defaults are
  package-dependent.
* `dtmin`: Minimum dt for adaptive timestepping. Defaults are
  package-dependent.
* `force_dtmin`: Declares whether to continue, forcing the minimum `dt` usage.
  Default is `false`, which has the solver throw a warning and exit early when
  encountering the minimum `dt`. Setting this true allows the solver to continue,
  never letting `dt` go below `dtmin` (and ignoring error tolerances in those
  cases). Note that `true` is not compatible with most interop packages.

#### Fixed Stepsize Usage

Note that if a method does not have adaptivity, the following rules apply:

* If `dt` is set, then the algorithm will step with size `dt` each iteration.
* If `tstops` and `dt` are both set, then the algorithm will step with either a
  size `dt`, or use a smaller step to hit the `tstops` point.
* If `tstops` is set without `dt`, then the algorithm will step directly to
  each value in `tstops`
* If neither `dt` nor `tstops` are set, the solver will throw an error.

#### [Advanced Adaptive Stepsize Control](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/)

These arguments control more advanced parts of the internals of adaptive timestepping
and are mostly used to make it more efficient on specific problems. For detailed
explanations of the timestepping algorithms, see the
[timestepping descriptions](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#timestepping)

* `internalnorm`: The norm function `internalnorm(u,t)` which error estimates
  are calculated. Required are two dispatches: one dispatch for the state variable
  and the other on the elements of the state variable (scalar norm).
  Defaults are package-dependent.
* `controller`: Possible examples are [`IController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.IController),
  [`PIController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.PIController),
  [`PIDController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.PIDController),
  [`PredictiveController`](https://docs.sciml.ai/DiffEqDocs/stable/extras/timestepping/#OrdinaryDiffEq.PredictiveController).
  Default is algorithm-dependent.
* `gamma`: The risk-factor γ in the q equation for adaptive timestepping
  of the controllers using it.
  Default is algorithm-dependent.
* `beta1`: The Lund stabilization α parameter.
  Default is algorithm-dependent.
* `beta2`: The Lund stabilization β parameter.
  Default is algorithm-dependent.
* `qmax`: Defines the maximum value possible for the adaptive q.
  Default is algorithm-dependent.
* `qmin`: Defines the minimum value possible for the adaptive q.
  Default is algorithm-dependent.
* `qsteady_min`: Defines the minimum for the range around 1 where the timestep
  is held constant. Default is algorithm-dependent.
* `qsteady_max`: Defines the maximum for the range around 1 where the timestep
  is held constant. Default is algorithm-dependent.
* `qoldinit`: The initial `qold` in stabilization stepping.
  Default is algorithm-dependent.
* `failfactor`: The amount to decrease the timestep by if the Newton iterations
  of an implicit method fail. Default is 2.

### Memory Optimizations

* `calck`: Turns on and off the internal ability for intermediate
  interpolations (also known as intermediate density). Not the same as `dense`, which is post-solution interpolation.
  This defaults to `dense || !isempty(saveat) ||  "no custom callback is given"`.
  This can be used to turn off interpolations
  (to save memory) if one isn't using interpolations when a custom callback is
  used. Another case where this may be used is to turn on interpolations for
  usage in the integrator interface even when interpolations are used nowhere else.
  Note that this is only required if the algorithm doesn't have
  a free or lazy interpolation (`DP8()`). If `calck = false`, `saveat` cannot be used.
  The rare keyword `calck` can be useful in event handling.
* `alias`: an `AbstractAliasSpecifier` object that holds fields specifying which variables to alias
  when solving. For example, to tell an ODE solver to alias the `u0` array, you can use an `ODEAliases` object, 
  and the `alias_u0` keyword argument, e.g. `solve(prob,alias = ODEAliases(alias_u0 = true))`. 
  For more information on what can be aliased for each problem type, see the documentation for the `AbstractAliasSpecifier`
  associated with that problem type. Set to `true` to alias every variable possible, or to `false` to disable aliasing.
  Defaults to an `AbstractAliasSpecifier` instance with `nothing` for all fields, which tells the solver to use the default behavior.

### Miscellaneous

* `maxiters`: Maximum number of iterations before stopping. Defaults to 1e5.
* `callback`: Specifies a callback. Defaults to a callback function which
  performs the saving routine. For more information, see the
  [Event Handling and Callback Functions manual page](https://docs.sciml.ai/DiffEqCallbacks/stable/).
* `initializealg`: The initialization algorithm for DAEs and ODEs with constraints.
  Available options include:
  - `DefaultInit()` (default): Automatically chooses the best initialization algorithm
  - `CheckInit()`: Only checks that initial conditions are consistent, errors if not
  - `NoInit()`: Skip initialization completely (for when you know conditions are consistent)
  - `OverrideInit()`: Use problem's initialization_data (typically from ModelingToolkit)
  - `BrownBasicInit()`: Brown's basic initialization algorithm for index-1 DAEs
  - `ShampineCollocationInit()`: Shampine's collocation initialization for general DAEs
  See the [DAE initialization documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/dae_initialization/)
  for more details.
* `isoutofdomain`: Specifies a function `isoutofdomain(u,p,t)` where, when it
  returns true, it will reject the timestep. Disabled by default.
* `unstable_check`: Specifies a function `unstable_check(dt,u,p,t)` where, when
  it returns true, it will cause the solver to exit and throw a warning. Defaults
  to `any(isnan,u)`, i.e. checking if any value is a NaN.
* `verbose`: Toggles whether warnings are thrown when the solver exits early.
  Defaults to true.
* `merge_callbacks`: Toggles whether to merge `prob.callback` with the `solve` keyword
  argument `callback`. Defaults to `true`.
* `wrap`: Toggles whether to wrap the solution if `prob.problem_type` has a preferred
  alternate wrapper type for the solution. Useful when speed, but not shape of solution
  is important. Defaults to `Val(true)`. `Val(false)` will cancel wrapping the solution.
* `u0`: The initial condition, overrides the one defined in the problem struct.
  Defaults to `nothing` (no override, use the `u0` defined in `prob`).
* `p`: The parameters, overrides the one defined in the problem struct.
  Defaults to `nothing` (no override, use the `p` defined in `prob`).

### Progress Monitoring

These arguments control the usage of the progressbar in ProgressLogging.jl compatible environments.
For information on setting up progress bars in VS Code and other environments, see the
[progress bar documentation](https://docs.sciml.ai/DiffEqDocs/stable/features/progress_bar/#Using-Progress-Bars-with-VS-Code).

* `progress`: Turns on/off the Juno progressbar. Default is false.
* `progress_steps`: Numbers of steps between updates of the progress bar.
  Default is 1000.
* `progress_name`: Controls the name of the progressbar. Default is the name
  of the problem type.
* `progress_message`: Controls the message with the progressbar. Defaults to
  showing `dt`, `t`, the maximum of `u`.
* `progress_id`: Controls the ID of the progress log message to distinguish simultaneous simulations.

### Error Calculations

If you are using the test problems (ex: `ODETestProblem`), then the following
options control the errors which are calculated:

* `timeseries_errors`: Turns on and off the calculation of errors at the steps
  which were taken, such as the `l2` error. Default is true.
* `dense_errors`: Turns on and off the calculation of errors at the steps which
  require dense output and calculate the error at 100 evenly-spaced points
  throughout `tspan`. An example is the `L2` error. Default is false.

### Sensitivity Algorithms (`sensealg`)

`sensealg` is used for choosing the way the automatic differentiation is performed.
For more information, see the documentation for SciMLSensitivity:
https://docs.sciml.ai/SciMLSensitivity/stable/

## Examples

The following lines are examples of how one could use the configuration of
`solve()`. For these examples a 3-dimensional ODE problem is assumed, however
the extension to other types is straightforward.

1. `solve(prob, AlgorithmName())` : The "default" setting, with a user-specified
   algorithm (given by `AlgorithmName()`). All parameters get their default values.
   This means that the solution is saved at the steps the Algorithm stops internally
   and dense output is enabled if the chosen algorithm allows for it.

   All other integration parameters (e.g. stepsize) are chosen automatically.
2. `solve(prob, saveat = 0.01, abstol = 1e-9, reltol = 1e-9)` : Standard setting
   for accurate output at specified (and equidistant) time intervals, used for
   e.g. Fourier Transform. The solution is given every 0.01 time units,
   starting from `tspan[1]`. The solver used is `Tsit5()` since no keyword
   `alg_hints` is given.

3. `solve(prob, maxiters = 1e7, progress = true, save_idxs = [1])` : Using longer
   maximum number of solver iterations can be useful when a given `tspan` is very
   long. This example only saves the first of the variables of the system, either
   to save size or because the user does not care about the others. Finally, with
   `progress = true` you are enabling the progress bar.
"""
function solve(
        prob::AbstractDEProblem, args...; sensealg = nothing,
        u0 = nothing, p = nothing, wrap = Val(true), kwargs...
    )
    if sensealg === nothing && haskey(prob.kwargs, :sensealg)
        sensealg = prob.kwargs[:sensealg]
    end

    u0 = u0 !== nothing ? u0 : prob.u0
    p = p !== nothing ? p : prob.p

    return if wrap isa Val{true}
        wrap_sol(
            solve_up(
                prob, sensealg, u0, p, args...;
                originator = SciMLBase.set_mooncakeoriginator_if_mooncake(SciMLBase.ChainRulesOriginator()),
                kwargs...
            )
        )
    else
        solve_up(
            prob, sensealg, u0, p, args...;
            originator = SciMLBase.set_mooncakeoriginator_if_mooncake(SciMLBase.ChainRulesOriginator()),
            kwargs...
        )
    end
end

function solve_up(
        prob::AbstractDEProblem, sensealg, u0, p,
        args...; originator = SciMLBase.ChainRulesOriginator(),
        kwargs...
    )
    alg = extract_alg(args, kwargs, has_kwargs(prob) ? prob.kwargs : kwargs)
    return if isnothing(alg) || !(alg isa AbstractDEAlgorithm) # Default algorithm handling
        _prob = get_concrete_problem(
            prob, !(prob isa DiscreteProblem); alg = alg, u0 = u0,
            p = p, kwargs...
        )
        solve_call(_prob, args...; kwargs...)
    else
        tstops = get(kwargs, :tstops, nothing)
        if tstops === nothing && has_kwargs(prob)
            tstops = get(prob.kwargs, :tstops, nothing)
        end
        if !(tstops isa Union{Nothing, AbstractArray, Tuple, Real}) &&
                !SciMLBase.allows_late_binding_tstops(alg)
            throw(LateBindingTstopsNotSupportedError())
        end
        _prob = get_concrete_problem(prob, isadaptive(alg); alg = alg, u0 = u0, p = p, kwargs...)
        _alg = prepare_alg(alg, _prob.u0, _prob.p, _prob)
        check_prob_alg_pairing(_prob, alg) # use alg for improved inference
        if length(args) > 1
            solve_call(_prob, _alg, Base.tail(args)...; kwargs...)
        else
            solve_call(_prob, _alg; kwargs...)
        end
    end
end

function solve(prob::AbstractNoiseProblem, args...; kwargs...)
    return __solve(prob, args...; kwargs...)
end

function solve(prob::AbstractJumpProblem, args...; kwargs...)
    return __solve(prob, args...; kwargs...)
end

function get_concrete_problem(prob::AbstractJumpProblem, isadapt; kwargs...)
    return get_updated_symbolic_problem(SciMLBase.get_root_indp(prob), prob; kwargs...)
end

function get_concrete_problem(prob::AbstractEnsembleProblem, isadapt; kwargs...)
    return prob
end

function solve(
        prob::PDEProblem, alg::AbstractDEAlgorithm, args...;
        kwargs...
    )
    return solve(prob.prob, alg, args...; kwargs...)
end

function init(
        prob::PDEProblem, alg::AbstractDEAlgorithm, args...;
        kwargs...
    )
    return init(prob.prob, alg, args...; kwargs...)
end

function get_concrete_problem(prob, isadapt; alg = nothing, kwargs...)
    oldprob = prob
    prob = get_updated_symbolic_problem(SciMLBase.get_root_indp(prob), prob; kwargs...)
    if prob !== oldprob
        kwargs = (; kwargs..., u0 = SII.state_values(prob), p = SII.parameter_values(prob))
    end
    p = get_concrete_p(prob, kwargs)
    tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
    u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
    u0_promote = promote_u0(u0, p, tspan[1])
    tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)
    f_promote = promote_f(
        prob.f, Val(SciMLBase.specialization(prob.f)), u0_promote, p,
        tspan_promote[1], Val(_uses_forwarddiff(alg)),
        _forwarddiff_chunksize(alg)
    )
    if isconcreteu0(prob, tspan[1], kwargs) && prob.u0 === u0 &&
            typeof(u0_promote) === typeof(prob.u0) &&
            prob.tspan == tspan && typeof(prob.tspan) === typeof(tspan_promote) &&
            p === prob.p && f_promote === prob.f
        return prob
    else
        return remake(prob; f = f_promote, u0 = u0_promote, p = p, tspan = tspan_promote)
    end
end

function get_concrete_problem(prob::DAEProblem, isadapt; alg = nothing, kwargs...)
    oldprob = prob
    prob = get_updated_symbolic_problem(SciMLBase.get_root_indp(prob), prob; kwargs...)
    if prob !== oldprob
        kwargs = (; kwargs..., u0 = SII.state_values(prob), p = SII.parameter_values(prob))
    end
    p = get_concrete_p(prob, kwargs)
    tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
    u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)
    du0 = get_concrete_du0(prob, isadapt, tspan[1], kwargs)

    u0_promote = promote_u0(u0, p, tspan[1])
    du0_promote = promote_u0(du0, p, tspan[1])
    tspan_promote = promote_tspan(u0_promote, p, tspan, prob, kwargs)

    f_promote = promote_f(
        prob.f, Val(SciMLBase.specialization(prob.f)), u0_promote, p,
        tspan_promote[1], Val(_uses_forwarddiff(alg)),
        _forwarddiff_chunksize(alg)
    )
    if isconcreteu0(prob, tspan[1], kwargs) && typeof(u0_promote) === typeof(prob.u0) &&
            isconcretedu0(prob, tspan[1], kwargs) && typeof(du0_promote) === typeof(prob.du0) &&
            prob.tspan == tspan && typeof(prob.tspan) === typeof(tspan_promote) &&
            p === prob.p && f_promote === prob.f
        return prob
    else
        return remake(
            prob; f = f_promote, du0 = du0_promote, u0 = u0_promote, p = p,
            tspan = tspan_promote
        )
    end
end

function get_concrete_problem(prob::DDEProblem, isadapt; kwargs...)
    oldprob = prob
    prob = get_updated_symbolic_problem(SciMLBase.get_root_indp(prob), prob; kwargs...)
    if prob !== oldprob
        kwargs = (; kwargs..., u0 = SII.state_values(prob), p = SII.parameter_values(prob))
    end
    p = get_concrete_p(prob, kwargs)
    tspan = get_concrete_tspan(prob, isadapt, kwargs, p)
    u0 = get_concrete_u0(prob, isadapt, tspan[1], kwargs)

    if prob.constant_lags isa Function
        constant_lags = prob.constant_lags(p)
    else
        constant_lags = prob.constant_lags
    end

    u0 = promote_u0(u0, p, tspan[1])
    tspan = promote_tspan(u0, p, tspan, prob, kwargs)

    return remake(prob; u0 = u0, tspan = tspan, p = p, constant_lags = constant_lags)
end

# Most are extensions
promote_tspan(u0, p, tspan, prob, kwargs) = _promote_tspan(tspan, kwargs)
function _promote_tspan(tspan, kwargs)
    if (dt = get(kwargs, :dt, nothing)) !== nothing
        tspan1, tspan2, _ = promote(tspan..., dt)
        return (tspan1, tspan2)
    else
        return tspan
    end
end

# Helper to get the effective ForwardDiff chunk size from the algorithm.
# Returns Val{CS}(). Defaults to Val(1) when algorithm is not known or unspecified.
_forwarddiff_chunksize(::Nothing) = Val(1)
_forwarddiff_chunksize(alg) = _resolve_chunksize(SciMLBase.forwarddiff_chunksize(alg))
_resolve_chunksize(::Val{0}) = Val(1)
_resolve_chunksize(v::Val) = v

# Helper to determine if we need ForwardDiff-aware function wrapping.
# Default to true (full wrapping) when algorithm is not known.
_uses_forwarddiff(::Nothing) = true
function _uses_forwarddiff(alg)
    fd = SciMLBase.forwarddiffs_model(alg)
    # forwarddiffs_model returns Bool for OrdinaryDiffEq algorithms, but may return
    # an AD type (e.g. AutoFiniteDiff()) for StochasticDiffEq algorithms.
    # Only use simple path when it explicitly returns false.
    fd !== false && return true
    # For composite algorithms that contain sub-algorithms (e.g. AutoTsit5(Rosenbrock23())),
    # check if any sub-algorithm uses ForwardDiff.
    if hasfield(typeof(alg), :algs)
        return any(_uses_forwarddiff, alg.algs)
    end
    return false
end

# Full path for algorithms that use ForwardDiff internally (e.g. Rosenbrock).
# These algorithms precompile AFTER the ForwardDiff extension loads, so
# backedges to hasdualpromote/wrapfun_iip don't cause invalidation issues.
function promote_f(
        f::F, ::Val{specialize}, u0, p, t, ::Val{true},
        ::Val{CS} = Val(1)
    ) where {F, specialize, CS}
    uElType = u0 === nothing ? Float64 : eltype(u0)
    if isdefined(f, :jac_prototype) && f.jac_prototype isa AbstractArray
        f = @set f.jac_prototype = similar(f.jac_prototype, uElType)
    end

    return f = if f isa ODEFunction && isinplace(f) && !(f.f isa AbstractSciMLOperator) &&
            # Opt-out SubArrays since they would create type mismatches with the integrator's internal Arrays
            !(u0 isa SubArray) &&
            (
            (
                specialize === SciMLBase.AutoSpecialize && eltype(u0) !== Any &&
                    RecursiveArrayTools.recursive_unitless_eltype(u0) === eltype(u0) &&
                    one(t) === oneunit(t) &&
                    hasdualpromote(u0, t)
            ) ||
                (
                specialize === SciMLBase.FunctionWrapperSpecialize &&
                    !(f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
            )
        )
        # Wrap tgrad if present, so its type is also erased.
        # tgrad!(dT, u, p, t) -> Nothing has the same shape as the RHS.
        if f.tgrad !== nothing && !(f.tgrad isa FunctionWrappersWrappers.FunctionWrappersWrapper)
            f = @set f.tgrad = wrapfun_jac_iip(f.tgrad, (u0, u0, p, t))
        end
        # Wrap the Jacobian if present, so its type is also erased.
        # Include both dense and sparse matrix signatures when the function
        # has a sparsity pattern, since the solver may use either depending on
        # the autodiff configuration (AutoSparse creates sparse J from sparsity).
        if f.jac !== nothing && !(f.jac isa FunctionWrappersWrappers.FunctionWrappersWrapper)
            if f.jac_prototype !== nothing
                J_T = Base.promote_op(similar, typeof(f.jac_prototype), Type{uElType})
                sig = Tuple{J_T, typeof(u0), typeof(p), typeof(t)}
                f = @set f.jac = FunctionWrappersWrappers.FunctionWrappersWrapper(
                    Void(f.jac), (sig,), (Nothing,))
            elseif isdefined(f, :sparsity) && f.sparsity isa AbstractMatrix &&
                    !(f.sparsity isa Matrix)
                # The sparsity pattern is a non-dense matrix (e.g. SparseMatrixCSC).
                # The solver may call the Jacobian with either a dense or sparse matrix
                # depending on the autodiff config, so wrap for both signatures.
                dense_sig = Tuple{Matrix{uElType}, typeof(u0), typeof(p), typeof(t)}
                sparse_J_T = Base.promote_op(similar, typeof(f.sparsity), Type{uElType})
                sparse_sig = Tuple{sparse_J_T, typeof(u0), typeof(p), typeof(t)}
                f = @set f.jac = FunctionWrappersWrappers.FunctionWrappersWrapper(
                    Void(f.jac),
                    (dense_sig, sparse_sig),
                    (Nothing, Nothing)
                )
            else
                sig = Tuple{Matrix{uElType}, typeof(u0), typeof(p), typeof(t)}
                f = @set f.jac = FunctionWrappersWrappers.FunctionWrappersWrapper(
                    Void(f.jac), (sig,), (Nothing,))
            end
        end
        return unwrapped_f(f, wrapfun_iip(f.f, (u0, u0, p, t), Val(CS)))
    else
        return f
    end
end

# Simple path for algorithms that do NOT use ForwardDiff internally (e.g. Tsit5, Verner).
# Avoids calling hasdualpromote/wrapfun_iip which have extension overrides in
# DiffEqBaseForwardDiffExt that would create invalidating method table backedges.
# Uses a simple single-signature FunctionWrapper instead.
function promote_f(
        f::F, ::Val{specialize}, u0, p, t, ::Val{false},
        ::Val{CS} = Val(1)
    ) where {F, specialize, CS}
    uElType = u0 === nothing ? Float64 : eltype(u0)
    if isdefined(f, :jac_prototype) && f.jac_prototype isa AbstractArray
        f = @set f.jac_prototype = similar(f.jac_prototype, uElType)
    end

    return f = if f isa ODEFunction && isinplace(f) && !(f.f isa AbstractSciMLOperator) &&
            f.mass_matrix isa UniformScaling &&
            f.jac === nothing &&
            !(u0 isa SubArray) &&
            (
            (
                specialize === SciMLBase.AutoSpecialize && eltype(u0) !== Any &&
                    RecursiveArrayTools.recursive_unitless_eltype(u0) === eltype(u0) &&
                    one(t) === oneunit(t)
            ) ||
                (
                specialize === SciMLBase.FunctionWrapperSpecialize &&
                    !(f.f isa FunctionWrappersWrappers.FunctionWrappersWrapper)
            )
        )
        return unwrapped_f(
            f,
            FunctionWrappersWrappers.FunctionWrappersWrapper(
                Void(f.f), (typeof((u0, u0, p, t)),), (Nothing,)
            )
        )
    else
        return f
    end
end

hasdualpromote(u0, t) = true

function promote_f(
        f::SplitFunction, ::Val{specialize}, u0, p, t, ::Val{true},
        ::Val{CS} = Val(1)
    ) where {specialize, CS}
    return if isnothing(f._func_cache)
        f
    else
        # Copy the cache to ensure it's properly initialized
        remake(f, _func_cache = copy(f._func_cache))
    end
end
function promote_f(
        f::SplitFunction, ::Val{specialize}, u0, p, t, ::Val{false},
        ::Val{CS} = Val(1)
    ) where {specialize, CS}
    return if isnothing(f._func_cache)
        f
    else
        remake(f, _func_cache = copy(f._func_cache))
    end
end
prepare_alg(alg, u0, p, f) = alg

function get_concrete_tspan(prob, isadapt, kwargs, p)
    if prob.tspan isa Function
        tspan = prob.tspan(p)
    elseif haskey(kwargs, :tspan)
        tspan = kwargs[:tspan]
    elseif prob.tspan === (nothing, nothing)
        throw(NoTspanError())
    else
        tspan = prob.tspan
    end

    isadapt && eltype(tspan) <: Integer && (tspan = float.(tspan))

    any(isnan, tspan) && throw(NaNTspanError())

    return tspan
end

function __solve(
        prob::AbstractDEProblem, args...; default_set = false, second_time = false,
        kwargs...
    )
    return if second_time
        throw(NoDefaultAlgorithmError())
    elseif length(args) > 0 && !(
            first(args) isa
                Union{Nothing, AbstractDEAlgorithm, AbstractNonlinearAlgorithm}
        )
        throw(NonSolverError())
    else
        __solve(prob, nothing, args...; default_set = false, second_time = true, kwargs...)
    end
end

function __init(
        prob::AbstractDEProblem, args...; default_set = false, second_time = false,
        kwargs...
    )
    return if second_time
        throw(NoDefaultAlgorithmError())
    elseif length(args) > 0 && !(
            first(args) isa
                Union{Nothing, AbstractDEAlgorithm, AbstractNonlinearAlgorithm}
        )
        throw(NonSolverError())
    else
        __init(prob, nothing, args...; default_set = false, second_time = true, kwargs...)
    end
end

function check_prob_alg_pairing(prob, alg)
    if prob isa ODEProblem && !(alg isa AbstractODEAlgorithm) ||
            prob isa SDEProblem && !(alg isa AbstractSDEAlgorithm) ||
            prob isa SDDEProblem && !(alg isa AbstractSDEAlgorithm) ||
            prob isa DDEProblem && !(alg isa AbstractDDEAlgorithm) ||
            prob isa DAEProblem && !(alg isa AbstractDAEAlgorithm) ||
            prob isa BVProblem && !(alg isa AbstractBVPAlgorithm) ||
            prob isa SteadyStateProblem && !(alg isa AbstractSteadyStateAlgorithm)
        throw(ProblemSolverPairingError(prob, alg))
    end

    if isdefined(prob, :u0) && eltypedual(prob.u0) &&
            !SciMLBase.isautodifferentiable(alg)
        throw(DirectAutodiffError())
    end

    if prob isa SDEProblem && prob.noise_rate_prototype !== nothing &&
            prob.noise !== nothing &&
            size(prob.noise_rate_prototype, 2) != length(prob.noise.W[1])
        throw(
            NoiseSizeIncompatibilityError(
                size(prob.noise_rate_prototype, 2),
                length(prob.noise.W[1])
            )
        )
    end

    # Complex number support comes before arbitrary number support for a more direct
    # error message.
    if !SciMLBase.allowscomplex(alg)
        if isdefined(prob, :u0) &&
                RecursiveArrayTools.recursive_unitless_eltype(prob.u0) <: Complex
            throw(ComplexSupportError(alg))
        end
    end

    if isdefined(prob, :tspan) && eltype(prob.tspan) <: Complex
        throw(ComplexTspanError())
    end

    # Check for concrete element type so that the non-concrete case throws a better error
    return if !SciMLBase.allows_arbitrary_number_types(alg)
        if isdefined(prob, :u0)
            uType = RecursiveArrayTools.recursive_unitless_eltype(prob.u0)
            u0_as_initial_guess = (prob isa BVProblem) && (uType <: Vector)
            if Base.isconcretetype(uType) &&
                    !(uType <: Union{Float32, Float64, ComplexF32, ComplexF64}) &&
                    !u0_as_initial_guess
                throw(
                    GenericNumberTypeError(
                        alg,
                        isdefined(prob, :u0) ? typeof(prob.u0) :
                            nothing,
                        isdefined(prob, :tspan) ? typeof(prob.tspan) :
                            nothing
                    )
                )
            end
        end

        if isdefined(prob, :tspan)
            tType = eltype(prob.tspan)
            if Base.isconcretetype(tType) &&
                    !(tType <: Union{Float32, Float64, ComplexF32, ComplexF64})
                throw(
                    GenericNumberTypeError(
                        alg,
                        isdefined(prob, :u0) ? typeof(prob.u0) :
                            nothing,
                        isdefined(prob, :tspan) ? typeof(prob.tspan) :
                            nothing
                    )
                )
            end
        end
    end
end

################### Differentiation

"""
Ignores all adjoint definitions (i.e. `sensealg`) and proceeds to do standard
AD through the `solve` functions. Generally only used internally for implementing
discrete sensitivity algorithms.
"""
struct SensitivityADPassThrough <: AbstractDEAlgorithm end

###
### Legacy Dispatches to be Non-Breaking
###

@deprecate concrete_solve(
    prob::AbstractDEProblem,
    alg::Union{AbstractDEAlgorithm, Nothing},
    u0 = prob.u0, p = prob.p, args...; kwargs...
) solve(
    prob, alg,
    args...;
    u0 = u0,
    p = p,
    kwargs...
)

function _solve_adjoint(
        prob, sensealg, u0, p, originator, args...; merge_callbacks = true,
        kwargs...
    )
    alg = extract_alg(args, kwargs, prob.kwargs)
    if isnothing(alg) || !(alg isa AbstractDEAlgorithm) # Default algorithm handling
        _prob = get_concrete_problem(
            prob, !(prob isa DiscreteProblem); alg = alg, u0 = u0,
            p = p, kwargs...
        )
    else
        _prob = get_concrete_problem(prob, isadaptive(alg); alg = alg, u0 = u0, p = p, kwargs...)
    end

    # Merge problem kwargs with passed kwargs
    kwargs = merge_problem_kwargs(_prob; merge_callbacks, kwargs...)

    return if length(args) > 1
        _concrete_solve_adjoint(
            _prob, alg, sensealg, u0, p, originator,
            Base.tail(args)...; kwargs...
        )
    else
        _concrete_solve_adjoint(_prob, alg, sensealg, u0, p, originator; kwargs...)
    end
end

function _solve_forward(
        prob, sensealg, u0, p, originator, args...; merge_callbacks = true,
        kwargs...
    )
    alg = extract_alg(args, kwargs, prob.kwargs)
    if isnothing(alg) || !(alg isa AbstractDEAlgorithm) # Default algorithm handling
        _prob = get_concrete_problem(
            prob, !(prob isa DiscreteProblem); alg = alg, u0 = u0,
            p = p, kwargs...
        )
    else
        _prob = get_concrete_problem(prob, isadaptive(alg); alg = alg, u0 = u0, p = p, kwargs...)
    end

    # Merge problem kwargs with passed kwargs
    kwargs = merge_problem_kwargs(_prob; merge_callbacks, kwargs...)

    return if length(args) > 1
        _concrete_solve_forward(
            _prob, alg, sensealg, u0, p, originator,
            Base.tail(args)...; kwargs...
        )
    else
        _concrete_solve_forward(_prob, alg, sensealg, u0, p, originator; kwargs...)
    end
end
