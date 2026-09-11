"""
    GlobalAdjoint(alg; gtol=nothing, adjoint_alg=alg, sensealg=nothing, samples=2,
                  rng=Random.default_rng(), maxiters=6, safety=0.8,
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

When solved with a `gtol`, the returned solution carries the endpoint estimate
in its `global_error` field (`SciMLBase.has_global_error` is `true` for
`GlobalAdjoint`); `sol.global_error[end]` is the estimated endpoint error
2-norm.

Omit `gtol` when using the algorithm only with [`adjoint_error_estimate`](@ref).
`maxiters` limits the number of forward-solve refinements, while `safety`
controls each local-tolerance reduction. The `adjoint_abstol` and
`adjoint_reltol` constructor keywords control the reverse-time adjoint solves.

The adjoint machinery lives in a package extension: SciMLSensitivity must be
loaded to solve with `GlobalAdjoint`. `sensealg` selects the adjoint sensitivity
algorithm; the default (`nothing`) resolves to `InterpolatingAdjoint` with a
ForwardDiff Jacobian (`autojacvec = false`). A ForwardDiff Jacobian is used
rather than a vector-Jacobian product because the estimator differentiates a
right-hand side that reads the forward solution's dense interpolant, which the
VJP backends mishandle (see SciML/SciMLSensitivity.jl#1649). `adjoint_alg`
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
    gtol::G
    options::V
end

function GlobalAdjoint(
        alg;
        adjoint_alg = alg,
        sensealg = nothing,
        samples = 2,
        rng = Random.default_rng(),
        gtol = nothing,
        maxiters = 6,
        safety = 0.8,
        adjoint_abstol = gtol === nothing ? 1.0e-10 : min(gtol / 100, 1.0e-10),
        adjoint_reltol = gtol === nothing ? 1.0e-8 : min(gtol / 100, 1.0e-8)
    )
    samples isa Integer || throw(ArgumentError("samples must be an integer"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    (gtol === nothing || _positive_finite_real(gtol)) ||
        throw(ArgumentError("gtol must be a positive finite real number"))
    maxiters isa Integer && maxiters > 0 ||
        throw(ArgumentError("maxiters must be a positive integer"))
    safety isa Real && isfinite(safety) && 0 < safety < 1 ||
        throw(ArgumentError("safety must be between zero and one"))
    _validate_tolerances(adjoint_abstol, adjoint_reltol, "adjoint")
    options = (; maxiters = Int(maxiters), safety, adjoint_abstol, adjoint_reltol)
    return GlobalAdjoint(alg, adjoint_alg, sensealg, Int(samples), rng, gtol, options)
end

# Extension hooks: GlobalDiffEqSciMLSensitivityExt adds methods for these when
# SciMLSensitivity is loaded.
function _adjoint_defect_projection end
function _default_adjoint_sensealg end

_adjoint_ext_loaded() = !isempty(methods(_default_adjoint_sensealg))

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
    return alg.sensealg === nothing ? _default_adjoint_sensealg() : alg.sensealg
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

function _adjoint_error_estimate(sol, alg, directions; adjoint_abstol, adjoint_reltol)
    sensealg = _resolve_sensealg(alg)
    projections = map(directions) do direction
        _adjoint_defect_projection(
            sol, sensealg, alg.adjoint_alg, direction;
            abstol = adjoint_abstol, reltol = adjoint_reltol
        )
    end
    T = eltype(sol.prob.u0)
    sample_count = length(directions)
    factor = _sphere_projection_expectation(sample_count, T) /
        _sphere_projection_expectation(length(sol.prob.u0), T)
    return factor * LinearAlgebra.norm(projections)
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
            final_sol = SciMLBase.solve(
                prob, alg.alg, args...;
                abstol = local_abstol, reltol = local_reltol, kwargs...
            )
            return @set final_sol.global_error =
                _endpoint_global_error(final_sol, last_estimate)
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
