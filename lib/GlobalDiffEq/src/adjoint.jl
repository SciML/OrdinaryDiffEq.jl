"""
    GlobalAdjoint(alg; gtol=nothing, adjoint_alg=alg, sensealg=nothing, samples=2,
                  rng=Random.default_rng(), maxiters=6, safety=0.8)

Wrap an ODE algorithm with adjoint-based global error estimation and control.
Set the requested absolute endpoint error with the `gtol` constructor keyword:

```julia
using SciMLSensitivity, QuadGK
solve(prob, GlobalAdjoint(Tsit5(); gtol = 1.0e-6))
```

The algorithm solves the ODE with local tolerances, estimates the 2-norm of the
endpoint error using `samples` orthogonal random projections, and tightens the
local tolerances until the estimate is at most `gtol`. The estimate is
probabilistic for systems with more than one state; the default of two samples
is accurate within a factor of ten with probability greater than 99% under the
small-sample estimate in Cao and Petzold (2004).

Omit `gtol` when using the algorithm only with [`adjoint_error_estimate`](@ref).
`maxiters` limits the number of forward-solve refinements, while `safety`
controls each local-tolerance reduction. The `adjoint_abstol`,
`adjoint_reltol`, `quadrature_abstol`, and `quadrature_reltol` constructor
keywords control the numerical adjoints and defect quadrature.

The adjoint machinery lives in a package extension: both SciMLSensitivity and
QuadGK must be loaded to solve with `GlobalAdjoint`. `sensealg` must be a
`SciMLSensitivity.QuadratureAdjoint`; the default (`nothing`) resolves to
`SciMLSensitivity.QuadratureAdjoint(autojacvec = true)`, which uses
SciMLSensitivity's ForwardDiff-based state Jacobian construction. Parameters
are presented to SciMLSensitivity as fixed data because this estimator does not
differentiate with respect to them. `adjoint_alg` selects the solver for the
reverse-time adjoint problems.

This implementation supports forward-time, standard-mass-matrix ODEs with
real vector states and no callbacks. The wrapped solver must provide a dense
first-derivative interpolation.
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
        adjoint_reltol = gtol === nothing ? 1.0e-8 : min(gtol / 100, 1.0e-8),
        quadrature_abstol = gtol === nothing ? 0.0 : zero(gtol),
        quadrature_reltol = 1.0e-6
    )
    samples isa Integer || throw(ArgumentError("samples must be an integer"))
    samples > 0 || throw(ArgumentError("samples must be positive"))
    sensealg === nothing || _is_quadrature_adjoint(sensealg) ||
        throw(ArgumentError(_ADJOINT_SENSEALG_MESSAGE))
    (gtol === nothing || _positive_finite_real(gtol)) ||
        throw(ArgumentError("gtol must be a positive finite real number"))
    maxiters isa Integer && maxiters > 0 ||
        throw(ArgumentError("maxiters must be a positive integer"))
    safety isa Real && isfinite(safety) && 0 < safety < 1 ||
        throw(ArgumentError("safety must be between zero and one"))
    _validate_tolerances(adjoint_abstol, adjoint_reltol, "adjoint")
    _validate_tolerances(quadrature_abstol, quadrature_reltol, "quadrature")
    options = (;
        maxiters = Int(maxiters), safety,
        adjoint_abstol, adjoint_reltol,
        quadrature_abstol, quadrature_reltol,
    )
    return GlobalAdjoint(alg, adjoint_alg, sensealg, Int(samples), rng, gtol, options)
end

# Extension hooks: GlobalDiffEqSciMLSensitivityExt adds methods for these when
# both SciMLSensitivity and QuadGK are loaded.
function _adjoint_solution end
function _defect_projection end
function _default_quadrature_sensealg end
_is_quadrature_adjoint(sensealg) = false

const _ADJOINT_SENSEALG_MESSAGE = "sensealg must be a SciMLSensitivity.QuadratureAdjoint, \
    and both SciMLSensitivity and QuadGK must be loaded for GlobalAdjoint \
    (run `using SciMLSensitivity, QuadGK`)"

_adjoint_ext_loaded() = !isempty(methods(_default_quadrature_sensealg))

function _require_adjoint_ext()
    _adjoint_ext_loaded() || throw(
        ArgumentError(
            "GlobalAdjoint requires the SciMLSensitivity extension: run " *
                "`using SciMLSensitivity, QuadGK` before solving"
        )
    )
    return nothing
end

function _resolve_sensealg(alg::GlobalAdjoint)
    _require_adjoint_ext()
    sensealg = alg.sensealg === nothing ? _default_quadrature_sensealg() : alg.sensealg
    _is_quadrature_adjoint(sensealg) || throw(ArgumentError(_ADJOINT_SENSEALG_MESSAGE))
    return sensealg
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

struct _FixedParameter{P, T}
    value::P
    tunables::Vector{T}
end

SciMLStructures.isscimlstructure(::_FixedParameter) = true
SciMLStructures.ismutablescimlstructure(::_FixedParameter) = false
SciMLStructures.hasportion(::SciMLStructures.AbstractPortion, ::_FixedParameter) = false
SciMLStructures.hasportion(::SciMLStructures.Tunable, ::_FixedParameter) = true

function SciMLStructures.canonicalize(::SciMLStructures.Tunable, p::_FixedParameter)
    repack = _ -> p
    return p.tunables, repack, true
end

function SciMLStructures.canonicalize(
        ::SciMLStructures.AbstractPortion, ::_FixedParameter
    )
    return nothing, nothing, nothing
end

function SciMLStructures.replace(
        ::SciMLStructures.Tunable, p::_FixedParameter, tunables
    )
    isempty(tunables) || throw(DimensionMismatch("fixed parameters have no tunable values"))
    return p
end

function _fixed_parameter_problem(prob)
    _validate_adjoint_problem(prob)
    fixed_parameter = _FixedParameter(prob.p, eltype(prob.u0)[])
    original_f = SciMLBase.unwrapped_f(prob.f)
    wrapped_f = if SciMLBase.isinplace(prob)
        let f = original_f
            (du, u, p, t) -> f(du, u, p.value, t)
        end
    else
        let f = original_f
            (u, p, t) -> f(u, p.value, t)
        end
    end
    full_specialized_f = SciMLBase.ODEFunction{
        SciMLBase.isinplace(prob), SciMLBase.FullSpecialize,
    }(wrapped_f)
    return SciMLBase.remake(prob; f = full_specialized_f, p = fixed_parameter)
end

function _validate_adjoint_problem(prob)
    prob.u0 isa AbstractVector{<:AbstractFloat} ||
        throw(ArgumentError("GlobalAdjoint requires a real floating-point vector state"))
    isempty(prob.u0) && throw(ArgumentError("GlobalAdjoint requires a nonempty state"))
    prob.tspan[1] < prob.tspan[2] ||
        throw(ArgumentError("GlobalAdjoint requires a forward-time ODEProblem"))
    prob.f.mass_matrix == LinearAlgebra.I ||
        throw(ArgumentError("GlobalAdjoint currently requires the standard mass matrix"))
    problem_kwargs = values(prob.kwargs)
    if haskey(problem_kwargs, :callback) && problem_kwargs.callback !== nothing
        throw(ArgumentError("GlobalAdjoint does not currently support callbacks"))
    end
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
        adjoint_abstol, adjoint_reltol, quadrature_abstol, quadrature_reltol
    )
    sensealg = _resolve_sensealg(alg)
    projections = map(directions) do direction
        adjoint_sol = _adjoint_solution(
            sol, sensealg, alg.adjoint_alg, direction;
            abstol = adjoint_abstol, reltol = adjoint_reltol
        )
        _defect_projection(
            sol, adjoint_sol;
            abstol = quadrature_abstol, reltol = quadrature_reltol
        )
    end
    T = eltype(sol.prob.u0)
    sample_count = length(directions)
    factor = _sphere_projection_expectation(sample_count, T) /
        _sphere_projection_expectation(length(sol.prob.u0), T)
    return factor * LinearAlgebra.norm(projections)
end

"""
    adjoint_error_estimate(prob, alg::GlobalAdjoint; abstol=1e-6, reltol=1e-3, kwargs...)

Solve an ODE problem and estimate the 2-norm of its global error at the final
time. The estimate evaluates the dense interpolation defect and projects it
through reverse-time adjoints constructed by SciMLSensitivity. Both
SciMLSensitivity and QuadGK must be loaded to enable the implementation.

`abstol` and `reltol` control the forward solve. Keyword arguments
`adjoint_abstol`, `adjoint_reltol`, `quadrature_abstol`, and `quadrature_reltol`
control the numerical adjoint solves and defect quadrature. The random
directions come from `alg.rng`; their number is `alg.samples`, capped at the
state dimension.
"""
function adjoint_error_estimate(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalAdjoint, args...;
        abstol = 1.0e-6,
        reltol = 1.0e-3,
        adjoint_abstol = alg.options.adjoint_abstol,
        adjoint_reltol = alg.options.adjoint_reltol,
        quadrature_abstol = alg.options.quadrature_abstol,
        quadrature_reltol = alg.options.quadrature_reltol,
        kwargs...
    )
    _require_adjoint_ext()
    haskey(kwargs, :callback) &&
        throw(ArgumentError("adjoint_error_estimate does not currently support callbacks"))
    estimation_prob = _fixed_parameter_problem(prob)
    solve_kwargs = merge(
        (; kwargs...),
        (;
            dense = true,
            save_everystep = true,
            save_start = true,
            save_end = true,
            saveat = (),
            save_idxs = nothing,
        )
    )
    sol = SciMLBase.solve(
        estimation_prob, alg.alg, args...; abstol, reltol, solve_kwargs...
    )
    _validate_adjoint_solution(sol)
    directions = _orthogonal_directions(prob.u0, alg.samples, alg.rng)
    return _adjoint_error_estimate(
        sol, alg, directions;
        adjoint_abstol,
        adjoint_reltol,
        quadrature_abstol,
        quadrature_reltol
    )
end

"""
    global_error_estimate(prob, alg::GlobalAdjoint; kwargs...)

Alias for [`adjoint_error_estimate`](@ref).
"""
function global_error_estimate(
        prob::SciMLBase.AbstractODEProblem, alg::GlobalAdjoint, args...; kwargs...
    )
    return adjoint_error_estimate(prob, alg, args...; kwargs...)
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
    haskey(kwargs, :callback) &&
        throw(ArgumentError("GlobalAdjoint does not currently support callbacks"))

    gtol = alg.gtol
    options = alg.options
    estimation_prob = _fixed_parameter_problem(prob)
    directions = _orthogonal_directions(prob.u0, alg.samples, alg.rng)
    local_abstol = abstol
    local_reltol = reltol
    trial_kwargs = merge(
        (; kwargs...),
        (;
            dense = true,
            save_everystep = true,
            save_start = true,
            save_end = true,
            saveat = (),
            save_idxs = nothing,
        )
    )
    last_estimate = oftype(float(gtol), Inf)

    for _ in 1:options.maxiters
        trial_sol = SciMLBase.solve(
            estimation_prob, alg.alg, args...;
            abstol = local_abstol, reltol = local_reltol, trial_kwargs...
        )
        _validate_adjoint_solution(trial_sol)
        last_estimate = _adjoint_error_estimate(
            trial_sol, alg, directions;
            adjoint_abstol = options.adjoint_abstol,
            adjoint_reltol = options.adjoint_reltol,
            quadrature_abstol = options.quadrature_abstol,
            quadrature_reltol = options.quadrature_reltol
        )
        isfinite(last_estimate) ||
            throw(ErrorException("the adjoint global error estimate is not finite"))

        if last_estimate <= gtol
            return SciMLBase.solve(
                prob, alg.alg, args...;
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
