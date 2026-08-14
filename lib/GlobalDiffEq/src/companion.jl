# Shared machinery for the defect-companion global error estimators
# (GlobalErrorTransport and GlobalDefectCorrection): both solve an auxiliary
# ODE for the global error driven by the defect of the dense numerical
# solution, d(t) = f(P(t)) - P'(t), where P is the solution interpolant.

# Tolerance validation helpers, shared with the other global-error estimators
# in this sublibrary. Kept local to the companion machinery so this file is
# self-contained; if the adjoint estimator lands too, these fold into one copy.
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

function _validate_estimation_problem(prob, name)
    prob.u0 isa AbstractVector{<:AbstractFloat} ||
        throw(ArgumentError("$name requires a real floating-point vector state"))
    isempty(prob.u0) && throw(ArgumentError("$name requires a nonempty state"))
    prob.tspan[1] < prob.tspan[2] ||
        throw(ArgumentError("$name requires a forward-time ODEProblem"))
    prob.f.mass_matrix == LinearAlgebra.I ||
        throw(ArgumentError("$name currently requires the standard mass matrix"))
    problem_kwargs = values(prob.kwargs)
    if haskey(problem_kwargs, :callback) && problem_kwargs.callback !== nothing
        throw(ArgumentError("$name does not currently support callbacks"))
    end
    return nothing
end

function _validate_estimation_solution(sol, name)
    SciMLBase.successful_retcode(sol) ||
        throw(ErrorException("the forward solve failed with retcode $(sol.retcode)"))
    length(sol.t) >= 2 ||
        throw(ArgumentError("the forward solution must save at least its two endpoints"))
    # The interpolating estimator differentiates the dense output (P' = sol(t, Val{1}))
    # to form the defect. A solver without genuine dense output leaves the trivial
    # Linear/Constant fallback interpolation, whose derivative does not approximate the
    # method's defect and would silently corrupt the estimate, so require real dense
    # output rather than trusting the caller.
    (sol.dense && !(sol.interp isa Union{
        SciMLBase.ConstantInterpolation, SciMLBase.LinearInterpolation,
    })) || throw(
        ArgumentError(
            "$name requires a solver with dense (continuous) output that provides a " *
                "first-derivative interpolation; the wrapped solver produced none. Use a " *
                "dense solver (e.g. Tsit5, Vern7) and do not disable `dense`."
        )
    )
    return nothing
end

"""
    GlobalErrorMode

Strategy a companion global-error wrapper ([`GlobalErrorTransport`](@ref),
[`GlobalDefectCorrection`](@ref)) uses to estimate the error. One of
[`InterpolatingMode`](@ref) or [`SimultaneousMode`](@ref).
"""
@enum GlobalErrorMode InterpolatingMode SimultaneousMode

@doc """
    InterpolatingMode

Run the forward solve with dense output and `save_everystep`, keep that dense
solution, and differentiate its interpolant to form the defect. Solver-agnostic,
but stores the whole dense solution (`O(n_steps)` memory) and runs a second
solve; requires a solver with a genuine dense first-derivative interpolation.
""" InterpolatingMode

@doc """
    SimultaneousMode

Co-integrate the error estimate in a single forward pass, using each step's
local dense output to inject the defect, without storing the whole dense
solution (`O(1)` extra memory, no second solve). See
SciML/OrdinaryDiffEq.jl#4491.
""" SimultaneousMode

# Only the interpolating implementation ships in this PR; the simultaneous
# implementation is the follow-up (SciML/OrdinaryDiffEq.jl#4491).
function _validate_gee_mode(mode::GlobalErrorMode)
    mode === InterpolatingMode || throw(
        ArgumentError(
            "mode must be InterpolatingMode, the only implemented mode; " *
                "SimultaneousMode is planned (SciML/OrdinaryDiffEq.jl#4491)"
        )
    )
    return nothing
end

# Out-of-place view of the (possibly in-place) user RHS, safe for dual numbers.
function _oop_rhs(prob)
    f = SciMLBase.unwrapped_f(prob.f)
    if SciMLBase.isinplace(prob)
        return let f = f
            function (u, p, t)
                du = similar(u)
                f(du, u, p, t)
                return du
            end
        end
    else
        return f
    end
end

const _DENSE_SOLVE_KWARGS = (;
    dense = true,
    save_everystep = true,
    save_start = true,
    save_end = true,
    saveat = (),
    save_idxs = nothing,
)

# Endpoint-error refinement loop shared by the companion estimators: solve,
# estimate the endpoint global error, and tighten the local tolerances until
# the estimate is at most gtol; then redo the solve with the user's original
# saving options.
function _refine_to_gtol(
        estimator, prob, inner_alg, gtol, options, args...;
        abstol, reltol, kwargs...
    )
    local_abstol = abstol
    local_reltol = reltol
    last_estimate = oftype(float(gtol), Inf)

    for _ in 1:options.maxiters
        last_estimate = estimator(local_abstol, local_reltol)
        isfinite(last_estimate) ||
            throw(ErrorException("the global error estimate is not finite"))

        if last_estimate <= gtol
            return SciMLBase.solve(
                prob, inner_alg, args...;
                abstol = local_abstol, reltol = local_reltol, kwargs...
            )
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

function _companion_options(gtol, maxiters, safety, companion_abstol, companion_reltol)
    maxiters isa Integer && maxiters > 0 ||
        throw(ArgumentError("maxiters must be a positive integer"))
    safety isa Real && isfinite(safety) && 0 < safety < 1 ||
        throw(ArgumentError("safety must be between zero and one"))
    (gtol === nothing || _positive_finite_real(gtol)) ||
        throw(ArgumentError("gtol must be a positive finite real number"))
    _validate_tolerances(companion_abstol, companion_reltol, "companion")
    return (;
        maxiters = Int(maxiters), safety,
        companion_abstol, companion_reltol,
    )
end

# Dense forward solve + companion error solve, returning the endpoint global
# error 2-norm estimate. rhs_for(sol) builds the companion RHS from the dense
# forward solution.
function _companion_error_estimate(
        rhs_for, name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
    haskey(kwargs, :callback) &&
        throw(ArgumentError("$name does not currently support callbacks"))
    _validate_estimation_problem(prob, name)
    solve_kwargs = merge((; kwargs...), _DENSE_SOLVE_KWARGS)
    sol = SciMLBase.solve(prob, inner_alg, args...; abstol, reltol, solve_kwargs...)
    _validate_estimation_solution(sol, name)
    endpoint_error = _companion_endpoint(
        rhs_for(sol), sol, companion_alg;
        abstol = companion_abstol, reltol = companion_reltol
    )
    return LinearAlgebra.norm(endpoint_error)
end

# Solve the companion error ODE ε' = rhs(ε, t), ε(t0) = 0, over the forward
# solution's time span and return the endpoint error estimate ε(T).
function _companion_endpoint(rhs, sol, companion_alg; abstol, reltol)
    ε0 = zero(sol.prob.u0)
    companion_prob = SciMLBase.ODEProblem{false}(rhs, ε0, sol.prob.tspan)
    companion_sol = SciMLBase.solve(
        companion_prob, companion_alg;
        abstol, reltol, dense = false, save_everystep = false,
        save_start = false, save_end = true
    )
    SciMLBase.successful_retcode(companion_sol) || throw(
        ErrorException(
            "the companion error solve failed with retcode $(companion_sol.retcode)"
        )
    )
    return companion_sol.u[end]
end
