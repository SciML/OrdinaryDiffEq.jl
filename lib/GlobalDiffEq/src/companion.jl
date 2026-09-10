# Shared machinery for the defect-companion global error estimator
# ([`GlobalErrorEstimation`](@ref)). After a dense forward solve producing the
# interpolant P(t), the global error ε(t) ≈ y(t) - P(t) is estimated by
# integrating a companion ODE driven by the defect d(t) = f(P(t)) - P'(t) of the
# dense output. The `equation` selects the companion ODE and the `mode` selects
# how it is integrated (see GlobalErrorEquation / GlobalErrorMode below).

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

# The estimator differentiates the dense output (P' = interp(t, Val{1})) to form
# the defect. A solver without genuine dense output leaves the trivial
# Linear/Constant fallback interpolation, whose derivative does not approximate
# the method's defect and would silently corrupt the estimate, so require real
# dense output rather than trusting the caller.
function _has_dense_derivative(sol)
    return sol.dense && !(sol.interp isa Union{
        SciMLBase.ConstantInterpolation, SciMLBase.LinearInterpolation,
    })
end

function _validate_estimation_solution(sol, name)
    SciMLBase.successful_retcode(sol) ||
        throw(ErrorException("the forward solve failed with retcode $(sol.retcode)"))
    length(sol.t) >= 2 ||
        throw(ArgumentError("the forward solution must save at least its two endpoints"))
    _has_dense_derivative(sol) || throw(
        ArgumentError(
            "$name requires a solver with dense (continuous) output that provides a " *
                "first-derivative interpolation; the wrapped solver produced none. Use a " *
                "dense solver (e.g. Tsit5, Vern7) and do not disable `dense`."
        )
    )
    return nothing
end

# Same dense-output requirement, checked on the integrator's solution before
# stepping in SimultaneousMode.
function _validate_streaming_interpolation(sol, name)
    _has_dense_derivative(sol) || throw(
        ArgumentError(
            "$name with SimultaneousMode requires a solver with dense (continuous) " *
                "output that provides a first-derivative interpolation; the wrapped " *
                "solver produced none. Use a dense solver (e.g. Tsit5, Vern7) and do " *
                "not disable `dense`."
        )
    )
    return nothing
end

"""
    GlobalErrorEquation

Which companion ODE [`GlobalErrorEstimation`](@ref) integrates for the global
error. One of [`DefectCorrection`](@ref) or [`ErrorTransport`](@ref).
"""
@enum GlobalErrorEquation DefectCorrection ErrorTransport

@doc """
    DefectCorrection

Integrate the nonlinear correction equation `ε' = f(P + ε, p, t) - P'(t)`,
`ε(t₀) = 0`, where `P` is the numerical solution's dense interpolant. Requires
no Jacobian, only extra evaluations of `f` on the perturbed argument `P + ε`.
The "solving for the correction" technique of Zadunaisky (1976) and Dormand,
Duckers and Prince (1984).
""" DefectCorrection

@doc """
    ErrorTransport

Integrate the linearized error-transport (first variational) equation
`ε' = J(P, p, t) ε + d(t)`, `ε(t₀) = 0`, where `d(t) = f(P, p, t) - P'(t)` is the
dense-output defect and `J = ∂f/∂u` is applied matrix-free as a Jacobian-vector
product (via `autodiff`). The error-transport approach of Shampine (1986) and
Berzins (1988). It is the first-order linearization of [`DefectCorrection`](@ref).
""" ErrorTransport

"""
    GlobalErrorMode

How [`GlobalErrorEstimation`](@ref) integrates the companion ODE. One of
[`InterpolatingMode`](@ref) or [`SimultaneousMode`](@ref).
"""
@enum GlobalErrorMode InterpolatingMode SimultaneousMode

@doc """
    InterpolatingMode

Run the forward solve with dense output and `save_everystep`, keep that dense
solution, and integrate the companion ODE against its interpolant in a second
solve. Solver-agnostic, but stores the whole dense solution (`O(n_steps)` memory)
and runs a second solve.
""" InterpolatingMode

@doc """
    SimultaneousMode

Co-integrate the error estimate in a single forward pass: drive the forward
integrator one accepted step at a time and, over each just-completed step,
advance the companion ODE using that step's live dense output as the interpolant.
Stores no dense solution (`O(1)` extra memory) and runs no second solve. Gives an
estimate closely agreeing with [`InterpolatingMode`](@ref).
""" SimultaneousMode

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

# The stored dense solution spans step nodes, where the defect jumps, so it is
# queried with right-continuity; a live integrator is restricted to the current
# step (continuity is both irrelevant and an unsupported keyword there), so it is
# queried plainly.
_interp_value(interp, t, right_continuity) =
    right_continuity ? interp(t, continuity = :right) : interp(t)
_interp_deriv(interp, t, right_continuity) =
    right_continuity ? interp(t, Val{1}, continuity = :right) : interp(t, Val{1})

# Nonlinear defect-correction companion RHS: ε' = f(P + ε) - P'.
function _defect_correction_rhs(prob, interp, right_continuity)
    foop = _oop_rhs(prob)
    return let interp = interp, foop = foop, p = prob.p, rc = right_continuity
        function (ε, _, t)
            u = _interp_value(interp, t, rc)
            du = _interp_deriv(interp, t, rc)
            return foop(u + ε, p, t) - du
        end
    end
end

# Linearized error-transport companion RHS: ε' = J(P) ε + (f(P) - P'), with the
# Jacobian applied matrix-free as a JVP through `autodiff`.
function _error_transport_rhs(prob, interp, right_continuity, autodiff)
    foop = _oop_rhs(prob)
    return let interp = interp, foop = foop, p = prob.p, rc = right_continuity,
            backend = autodiff

        function (ε, _, t)
            u = _interp_value(interp, t, rc)
            du = _interp_deriv(interp, t, rc)
            defect = foop(u, p, t) - du
            jv = only(
                DifferentiationInterface.pushforward(
                    x -> foop(x, p, t), backend, u, (ε,)
                )
            )
            return jv + defect
        end
    end
end

# Return `make_rhs(interp, prob, right_continuity)` selecting the companion ODE.
function _companion_rhs_builder(equation::GlobalErrorEquation, autodiff)
    if equation === DefectCorrection
        return (interp, prob, rc) -> _defect_correction_rhs(prob, interp, rc)
    else
        return (interp, prob, rc) -> _error_transport_rhs(prob, interp, rc, autodiff)
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

# Endpoint-error refinement loop: solve, estimate the endpoint global error, and
# tighten the local tolerances until the estimate is at most gtol; then redo the
# solve with the user's original saving options.
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

# InterpolatingMode: dense forward solve, then integrate the companion against
# the stored interpolant. `make_rhs(interp, prob, right_continuity)` builds the
# companion RHS. Returns the endpoint global error 2-norm estimate.
function _companion_error_estimate(
        make_rhs, name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
    haskey(kwargs, :callback) &&
        throw(ArgumentError("$name does not currently support callbacks"))
    _validate_estimation_problem(prob, name)
    solve_kwargs = merge((; kwargs...), _DENSE_SOLVE_KWARGS)
    sol = SciMLBase.solve(prob, inner_alg, args...; abstol, reltol, solve_kwargs...)
    _validate_estimation_solution(sol, name)
    endpoint_error = _companion_endpoint(
        make_rhs(sol, sol.prob, true), sol, companion_alg;
        abstol = companion_abstol, reltol = companion_reltol
    )
    return LinearAlgebra.norm(endpoint_error)
end

# SimultaneousMode: drive the forward integrator one accepted step at a time and
# advance the companion across each just-completed step [tprev, t] using the
# integrator's live interpolant, keeping no dense solution and running no second
# solve. Returns the endpoint global error 2-norm estimate.
function _companion_error_estimate_streaming(
        make_rhs, name, prob, inner_alg, companion_alg, args...;
        abstol, reltol, companion_abstol, companion_reltol, kwargs...
    )
    haskey(kwargs, :callback) &&
        throw(ArgumentError("$name does not currently support callbacks"))
    _validate_estimation_problem(prob, name)
    integrator = SciMLBase.init(
        prob, inner_alg, args...;
        abstol, reltol, dense = true, save_everystep = false,
        save_start = true, save_end = true, kwargs...
    )
    _validate_streaming_interpolation(integrator.sol, name)
    ε = zero(prob.u0)
    for _ in integrator
        tprev, t = integrator.tprev, integrator.t
        t > tprev || continue
        ε = _advance_companion(
            make_rhs(integrator, prob, false), ε, tprev, t, companion_alg;
            abstol = companion_abstol, reltol = companion_reltol
        )
    end
    SciMLBase.successful_retcode(integrator.sol) ||
        throw(ErrorException("the forward solve failed with retcode $(integrator.sol.retcode)"))
    return LinearAlgebra.norm(ε)
end

# Solve the companion ODE ε' = rhs(ε, t), ε(t0) = 0, over the forward solution's
# time span and return the endpoint error estimate ε(T).
function _companion_endpoint(rhs, sol, companion_alg; abstol, reltol)
    return _advance_companion(
        rhs, zero(sol.prob.u0), sol.prob.tspan[1], sol.prob.tspan[2],
        companion_alg; abstol, reltol
    )
end

# Advance the companion ODE ε' = rhs(ε, t) from ε(t0) = ε0 over [t0, t1] and
# return ε(t1).
function _advance_companion(rhs, ε0, t0, t1, companion_alg; abstol, reltol)
    companion_prob = SciMLBase.ODEProblem{false}(rhs, ε0, (t0, t1))
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
