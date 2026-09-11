"""
    GlobalAdjointScope

What [`GlobalAdjoint`](@ref) estimates and reports in `sol.global_error`. One of
[`EndpointError`](@ref) or [`TrajectoryError`](@ref).
"""
@enum GlobalAdjointScope EndpointError TrajectoryError

@doc """
    EndpointError

Estimate the global error 2-norm only at the final time: `sol.global_error[end]`
holds it and earlier entries are zero. The default; one adjoint solve per sample.
""" EndpointError

@doc """
    TrajectoryError

Estimate the global error 2-norm at every saved time: `sol.global_error[i]` is the
estimate at `sol.t[i]` (`sol.global_error[1] = 0`), filled along the whole
trajectory. Costs one reverse-adjoint solve per saved time per sample, and the
returned solution is dense and saved at every solver step.
""" TrajectoryError

"""
    GlobalAdjointControl

How [`GlobalAdjoint`](@ref) drives the solve to meet `gtol`. One of
[`ToleranceRefinement`](@ref) or [`StepGridRefinement`](@ref).
"""
@enum GlobalAdjointControl ToleranceRefinement StepGridRefinement

@doc """
    ToleranceRefinement

Tighten the solver's local `abstol`/`reltol` and re-solve adaptively until the
estimated endpoint error is at most `gtol`. The default.
""" ToleranceRefinement

@doc """
    StepGridRefinement

Tighten tolerances the same way to find a step grid that meets `gtol`, then return
a solution obtained by re-solving **non-adaptively on exactly that grid**
(`adaptive = false`, stepping on the adjoint-informed `dt`s). Yields a
reproducible fixed-step schedule rather than an adaptive solve.
""" StepGridRefinement

"""
    GlobalAdjoint(alg; gtol=nothing, adjoint_alg=alg, sensealg=nothing, samples=2,
                  rng=Random.default_rng(), scope=EndpointError,
                  control=ToleranceRefinement, maxiters=6, safety=0.8,
                  adjoint_abstol, adjoint_reltol)

Wrap an ODE algorithm with adjoint-based global error estimation and control.
Set the requested absolute endpoint error with the `gtol` constructor keyword:

```julia
using SciMLSensitivity
solve(prob, GlobalAdjoint(Tsit5(); gtol = 1.0e-6))
```

The algorithm solves the ODE with local tolerances, estimates the 2-norm of the
endpoint error using `samples` orthogonal random projections, and tightens the
local tolerances until the estimate is at most `gtol`. The estimate is
probabilistic for systems with more than one state; the default of two samples
is accurate within a factor of ten with probability greater than 99% under the
small-sample estimate in Cao and Petzold (2004).

The estimate is a dual-weighted residual: for a projection direction `v`, the
error component `⟨v, y(T) - P(T)⟩` equals `∫ λ(t)ᵀ (f(P(t)) - P'(t)) dt`, where
`P` is the dense numerical solution, `d(t) = f(P) - P'` its defect, and `λ` the
adjoint with terminal condition `λ(T) = v`. This is computed by the adjoint
system itself (via `SciMLSensitivity.adjoint_sensitivities`), not by a
hand-rolled quadrature: the defect is introduced as a scalar forcing whose
sensitivity is exactly that integral.

When solved with a `gtol`, the returned solution carries the estimate in its
`global_error` field (`SciMLBase.has_global_error` is `true` for `GlobalAdjoint`).
`scope` ([`GlobalAdjointScope`](@ref)) selects what is estimated:
[`EndpointError`](@ref) (default) fills only `sol.global_error[end]`, while
[`TrajectoryError`](@ref) fills `sol.global_error[i]` at every saved time.
`control` ([`GlobalAdjointControl`](@ref)) selects how `gtol` is met:
[`ToleranceRefinement`](@ref) (default) tightens tolerances and re-solves
adaptively, while [`StepGridRefinement`](@ref) re-solves non-adaptively on the
accepted step grid. Control always targets the endpoint error estimate.

Omit `gtol` when using the algorithm only with [`adjoint_error_estimate`](@ref).
`maxiters` limits the number of forward-solve refinements, while `safety`
controls each local-tolerance reduction. The `adjoint_abstol` and
`adjoint_reltol` constructor keywords control the reverse-time adjoint solves.

The adjoint machinery lives in a package extension: SciMLSensitivity must be
loaded to solve with `GlobalAdjoint`. `sensealg` selects the adjoint sensitivity
algorithm and is passed straight through to `SciMLSensitivity.adjoint_sensitivities`:
any adjoint method works, and the default (`nothing`) defers to that function's
own default. `adjoint_alg`
selects the solver for the reverse-time adjoint problems.

This implementation supports forward-time, standard-mass-matrix ODEs with real
vector states. The wrapped solver must provide a dense first-derivative
interpolation.

## References

  - Y. Cao and L. Petzold, A posteriori error estimation and global error
    control for ordinary differential equations by the adjoint method, SIAM
    Journal on Scientific Computing 26 (2004).
"""
struct GlobalAdjoint{A, AA, S, R, G, V} <: GlobalDiffEqAlgorithm
    alg::A
    adjoint_alg::AA
    sensealg::S
    samples::Int
    rng::R
    scope::GlobalAdjointScope
    control::GlobalAdjointControl
    gtol::G
    options::V
end

function GlobalAdjoint(
        alg;
        adjoint_alg = alg,
        sensealg = nothing,
        samples = 2,
        rng = Random.default_rng(),
        scope = EndpointError,
        control = ToleranceRefinement,
        gtol = nothing,
        maxiters = 6,
        safety = 0.8,
        adjoint_abstol = gtol === nothing ? 1.0e-10 : min(gtol / 100, 1.0e-10),
        adjoint_reltol = gtol === nothing ? 1.0e-8 : min(gtol / 100, 1.0e-8)
    )
    samples isa Integer || throw(ArgumentError("samples must be an integer"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    scope isa GlobalAdjointScope ||
        throw(ArgumentError("scope must be EndpointError or TrajectoryError"))
    control isa GlobalAdjointControl ||
        throw(ArgumentError("control must be ToleranceRefinement or StepGridRefinement"))
    (gtol === nothing || _positive_finite_real(gtol)) ||
        throw(ArgumentError("gtol must be a positive finite real number"))
    maxiters isa Integer && maxiters > 0 ||
        throw(ArgumentError("maxiters must be a positive integer"))
    safety isa Real && isfinite(safety) && 0 < safety < 1 ||
        throw(ArgumentError("safety must be between zero and one"))
    _validate_tolerances(adjoint_abstol, adjoint_reltol, "adjoint")
    options = (; maxiters = Int(maxiters), safety, adjoint_abstol, adjoint_reltol)
    return GlobalAdjoint(
        alg, adjoint_alg, sensealg, Int(samples), rng, scope, control, gtol, options
    )
end

# Extension hook: GlobalDiffEqSciMLSensitivityExt adds a method for this when
# SciMLSensitivity is loaded.
function _adjoint_defect_projection end

_adjoint_ext_loaded() = !isempty(methods(_adjoint_defect_projection))

function _require_adjoint_ext()
    _adjoint_ext_loaded() || throw(
        ArgumentError(
            "GlobalAdjoint requires the SciMLSensitivity extension: run " *
                "`using SciMLSensitivity` before solving"
        )
    )
    return nothing
end

function _resolve_sensealg(alg::GlobalAdjoint)
    _require_adjoint_ext()
    return alg.sensealg
end

_positive_finite_real(value) = value isa Real && isfinite(value) && value > 0

function _validate_tolerances(abstol, reltol, name)
    abstol isa Real && isfinite(abstol) && abstol >= 0 ||
        throw(ArgumentError("$(name)_abstol must be a nonnegative finite real number"))
    reltol isa Real && isfinite(reltol) && reltol >= 0 ||
        throw(ArgumentError("$(name)_reltol must be a nonnegative finite real number"))
    iszero(abstol) && iszero(reltol) &&
        throw(ArgumentError("$(name)_abstol and $(name)_reltol cannot both be zero"))
    return nothing
end

SciMLBase.allows_arbitrary_number_types(::GlobalAdjoint) = false
SciMLBase.allowscomplex(::GlobalAdjoint) = false
SciMLBase.isautodifferentiable(::GlobalAdjoint) = false
SciMLBase.has_global_error(::GlobalAdjoint) = true

function _validate_adjoint_problem(prob)
    prob.u0 isa AbstractVector{<:AbstractFloat} ||
        throw(ArgumentError("GlobalAdjoint requires a real floating-point vector state"))
    isempty(prob.u0) && throw(ArgumentError("GlobalAdjoint requires a nonempty state"))
    prob.tspan[1] < prob.tspan[2] ||
        throw(ArgumentError("GlobalAdjoint requires a forward-time ODEProblem"))
    prob.f.mass_matrix == LinearAlgebra.I ||
        throw(ArgumentError("GlobalAdjoint currently requires the standard mass matrix"))
    return nothing
end

function _validate_adjoint_solution(sol)
    SciMLBase.successful_retcode(sol) ||
        throw(ErrorException("the forward solve failed with retcode $(sol.retcode)"))
    _validate_adjoint_problem(sol.prob)
    length(sol.t) >= 2 ||
        throw(ArgumentError("the forward solution must save at least its two endpoints"))
    return nothing
end

function _orthogonal_directions(u0, requested_samples, rng)
    sample_count = min(requested_samples, length(u0))
    T = eltype(u0)
    directions = Vector{Vector{T}}()
    sizehint!(directions, sample_count)

    for _ in 1:sample_count
        accepted = false
        for _ in 1:10
            direction = T[Base.randn(rng) for _ in eachindex(u0)]
            for _ in 1:2, previous in directions
                direction .-= LinearAlgebra.dot(previous, direction) .* previous
            end
            direction_norm = LinearAlgebra.norm(direction)
            if direction_norm > sqrt(eps(T))
                push!(directions, direction ./ direction_norm)
                accepted = true
                break
            end
        end
        accepted || error("failed to generate orthogonal random directions")
    end

    return directions
end

function _sphere_projection_expectation(dimension, ::Type{T}) where {T}
    dimension > 0 || throw(ArgumentError("dimension must be positive"))
    if dimension == 1
        return one(T)
    end

    expectation = isodd(dimension) ? one(T) : T(2) / T(pi)
    first_numerator = isodd(dimension) ? 1 : 2
    for numerator in first_numerator:2:(dimension - 2)
        expectation *= T(numerator) / T(numerator + 1)
    end
    return expectation
end

function _adjoint_error_estimate(
        sol, alg, directions;
        adjoint_abstol, adjoint_reltol, terminal_time = sol.prob.tspan[2]
    )
    sensealg = _resolve_sensealg(alg)
    projections = map(directions) do direction
        _adjoint_defect_projection(
            sol, sensealg, alg.adjoint_alg, direction;
            abstol = adjoint_abstol, reltol = adjoint_reltol, terminal_time
        )
    end
    T = eltype(sol.prob.u0)
    sample_count = length(directions)
    factor = _sphere_projection_expectation(sample_count, T) /
        _sphere_projection_expectation(length(sol.prob.u0), T)
    return factor * LinearAlgebra.norm(projections)
end

# Per-time global error along the trajectory: the estimate at each saved time,
# obtained by placing the adjoint's discrete cost at that time. `sol.t[1]` (the
# initial time) carries no error.
function _adjoint_trajectory_errors(sol, alg, directions; adjoint_abstol, adjoint_reltol)
    errors = zeros(eltype(sol.prob.u0), length(sol.t))
    for j in 2:length(sol.t)
        errors[j] = _adjoint_error_estimate(
            sol, alg, directions;
            adjoint_abstol, adjoint_reltol, terminal_time = sol.t[j]
        )
    end
    return errors
end

# Re-solve non-adaptively on exactly `grid`, stepping on its `dt`s. Produces a
# dense solution saved at every grid point, matching the adaptive solve's steps.
function _solve_on_grid(prob, inner_alg, grid, args...; abstol, reltol)
    dts = diff(grid)
    integrator = SciMLBase.init(
        prob, inner_alg, args...;
        adaptive = false, dt = dts[1], dense = true,
        save_everystep = true, save_start = true, save_end = true,
        abstol, reltol
    )
    for dt in dts
        SciMLBase.set_proposed_dt!(integrator, dt)
        SciMLBase.step!(integrator)
    end
    sol = integrator.sol
    if !SciMLBase.successful_retcode(sol) && integrator.t >= grid[end]
        sol = @set sol.retcode = SciMLBase.ReturnCode.Success
    end
    return sol
end

const _DENSE_SOLVE_KWARGS = (;
    dense = true,
    save_everystep = true,
    save_start = true,
    save_end = true,
    saveat = (),
    save_idxs = nothing,
)

"""
    adjoint_error_estimate(prob, alg::GlobalAdjoint; abstol=1e-6, reltol=1e-3, kwargs...)

Solve an ODE problem and estimate the 2-norm of its global error at the final
time. The estimate evaluates the dense interpolation defect and projects it
through reverse-time adjoints constructed by SciMLSensitivity, which must be
loaded to enable the implementation.

`abstol` and `reltol` control the forward solve. `adjoint_abstol` and
`adjoint_reltol` control the reverse-time adjoint solves. The random directions
come from `alg.rng`; their number is `alg.samples`, capped at the state
dimension.
"""
function adjoint_error_estimate(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalAdjoint, args...;
        abstol = 1.0e-6,
        reltol = 1.0e-3,
        adjoint_abstol = alg.options.adjoint_abstol,
        adjoint_reltol = alg.options.adjoint_reltol,
        kwargs...
    )
    _require_adjoint_ext()
    _validate_adjoint_problem(prob)
    solve_kwargs = merge((; kwargs...), _DENSE_SOLVE_KWARGS)
    sol = SciMLBase.solve(prob, alg.alg, args...; abstol, reltol, solve_kwargs...)
    _validate_adjoint_solution(sol)
    directions = _orthogonal_directions(prob.u0, alg.samples, alg.rng)
    return _adjoint_error_estimate(sol, alg, directions; adjoint_abstol, adjoint_reltol)
end

# Endpoint-mode global error: a vector aligned with the saved times carrying the
# estimated endpoint 2-norm in the final slot (earlier times are not estimated in
# endpoint mode).
function _endpoint_global_error(sol, estimate)
    global_error = zeros(typeof(estimate), length(sol.t))
    global_error[end] = estimate
    return global_error
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalAdjoint, args...;
        abstol = something(alg.gtol, 1.0e-6),
        reltol = something(alg.gtol, 1.0e-3),
        kwargs...
    )
    alg.gtol === nothing &&
        throw(ArgumentError("GlobalAdjoint requires a positive `gtol` constructor keyword"))
    _require_adjoint_ext()
    _validate_tolerances(abstol, reltol, "local")
    _validate_adjoint_problem(prob)

    gtol = alg.gtol
    options = alg.options
    directions = _orthogonal_directions(prob.u0, alg.samples, alg.rng)
    local_abstol = abstol
    local_reltol = reltol
    trial_kwargs = merge((; kwargs...), _DENSE_SOLVE_KWARGS)
    last_estimate = oftype(float(gtol), Inf)

    for _ in 1:options.maxiters
        trial_sol = SciMLBase.solve(
            prob, alg.alg, args...;
            abstol = local_abstol, reltol = local_reltol, trial_kwargs...
        )
        _validate_adjoint_solution(trial_sol)
        last_estimate = _adjoint_error_estimate(
            trial_sol, alg, directions;
            adjoint_abstol = options.adjoint_abstol,
            adjoint_reltol = options.adjoint_reltol
        )
        isfinite(last_estimate) ||
            throw(ErrorException("the adjoint global error estimate is not finite"))

        if last_estimate <= gtol
            # The final solution reproduces the accepted grid non-adaptively
            # (StepGridRefinement) or re-solves adaptively with the user's saving
            # options (ToleranceRefinement). TrajectoryError needs the dense,
            # every-step solution to estimate along the whole path.
            final_sol = if alg.control === StepGridRefinement
                _solve_on_grid(
                    prob, alg.alg, trial_sol.t, args...;
                    abstol = local_abstol, reltol = local_reltol
                )
            elseif alg.scope === TrajectoryError
                SciMLBase.solve(
                    prob, alg.alg, args...;
                    abstol = local_abstol, reltol = local_reltol, trial_kwargs...
                )
            else
                SciMLBase.solve(
                    prob, alg.alg, args...;
                    abstol = local_abstol, reltol = local_reltol, kwargs...
                )
            end
            global_error = if alg.scope === TrajectoryError
                _adjoint_trajectory_errors(
                    final_sol, alg, directions;
                    adjoint_abstol = options.adjoint_abstol,
                    adjoint_reltol = options.adjoint_reltol
                )
            else
                _endpoint_global_error(final_sol, last_estimate)
            end
            return @set final_sol.global_error = global_error
        end

        tolerance_scale = min(0.5, options.safety * gtol / last_estimate)
        local_abstol *= tolerance_scale
        local_reltol *= tolerance_scale
        iszero(local_abstol) && iszero(local_reltol) &&
            throw(ErrorException("local tolerances underflowed during global error control"))
    end

    throw(
        ErrorException(
            "failed to meet gtol=$(gtol) after $(options.maxiters) iterations; " *
                "last estimated global error was $(last_estimate)"
        )
    )
end
