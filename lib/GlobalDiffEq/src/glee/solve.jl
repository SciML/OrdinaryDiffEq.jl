# The GLEE methods integrate a partitioned state (y, ε). Users pass a plain
# ODEProblem; __solve extends it to ArrayPartition form with this RHS wrapper,
# which applies the user's f to the solution partition. The ε partition has no
# ODE of its own (the GL method updates it directly), so its rate is zero
# wherever the wrapper is evaluated (initialization, dense output derivatives).
struct GLEEExtendedRHS{iip, F}
    f::F
end

function (rhs::GLEEExtendedRHS{true})(du, u, p, t)
    rhs.f(du.x[1], u.x[1], p, t)
    fill!(du.x[2], zero(eltype(du.x[2])))
    return nothing
end

function (rhs::GLEEExtendedRHS{false})(u, p, t)
    return RecursiveArrayTools.ArrayPartition(rhs.f(u.x[1], p, t), zero(u.x[1]))
end

function _glee_inner(f)
    rhs = f isa SciMLBase.AbstractODEFunction ? f.f : f
    rhs isa GLEEExtendedRHS || throw(
        ArgumentError(
            "GLEE methods integrate plain ODEProblems through their own " *
                "partitioned-state extension; do not pass a manually partitioned problem"
        )
    )
    return rhs.f
end

function _glee_extended_problem(prob)
    prob.u0 isa RecursiveArrayTools.ArrayPartition && throw(
        ArgumentError(
            "GLEE methods construct their own partitioned (y, ε) state; " *
                "pass the plain ODEProblem instead of an ArrayPartition state"
        )
    )
    prob.u0 isa AbstractArray ||
        throw(ArgumentError("GLEE methods require an array state"))
    prob.f.mass_matrix == LinearAlgebra.I ||
        throw(ArgumentError("GLEE methods require the standard mass matrix"))
    iip = SciMLBase.isinplace(prob)
    rhs = GLEEExtendedRHS{iip, typeof(prob.f)}(prob.f)
    extended_f = SciMLBase.ODEFunction{iip, SciMLBase.FullSpecialize}(rhs)
    u0 = RecursiveArrayTools.ArrayPartition(copy(prob.u0), zero(prob.u0))
    return SciMLBase.remake(prob; f = extended_f, u0)
end

_is_glee_extended(prob) = prob.f.f isa GLEEExtendedRHS

# init/solve on a plain ODEProblem transparently extend it to the partitioned
# (y, ε) state; the invoke dispatches into OrdinaryDiffEqCore's generic __init
# (its exact five-argument form, so the DiffEqBase default-algorithm catch-all
# cannot be selected).
function SciMLBase.__init(
        prob::SciMLBase.AbstractODEProblem, alg::AbstractGLEEAlgorithm,
        timeseries_init = (), ts_init = (), ks_init = ();
        kwargs...
    )
    extended_prob = _is_glee_extended(prob) ? prob : _glee_extended_problem(prob)
    return invoke(
        SciMLBase.__init,
        Tuple{
            SciMLBase.AbstractODEProblem,
            OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm,
            Any, Any, Any,
        },
        extended_prob, alg, timeseries_init, ts_init, ks_init; kwargs...
    )
end

# Wraps the interpolation of the extended (y, ε) solve and projects it onto the
# solution partition, so `sol(t)` returns the ordinary solution at full
# interpolation order while the global error estimate lives in
# `sol.global_error`.
struct GLEESolutionInterpolation{I} <: SciMLBase.AbstractDiffEqInterpolation
    inner::I
end

_project_y(u::RecursiveArrayTools.ArrayPartition) = u.x[1]
_project_y(u::AbstractVector) = map(_project_y, u)

function _project_interpolated(full, idxs)
    if full isa RecursiveArrayTools.AbstractDiffEqArray
        projected = RecursiveArrayTools.DiffEqArray(map(_project_y, full.u), full.t)
        return idxs === nothing ? projected : projected[idxs]
    end
    y = _project_y(full)
    return idxs === nothing ? y : y[idxs]
end

# The extended interpolation is always queried with `idxs = nothing` (whole
# partitioned state), then projected; user `idxs` index into the projected
# solution component, never into the raw partition layout.
function (interp::GLEESolutionInterpolation)(
        tvals, idxs, deriv, p, continuity::Symbol = :left
    )
    full = interp.inner(tvals, nothing, deriv, p, continuity)
    return _project_interpolated(full, idxs)
end

function (interp::GLEESolutionInterpolation)(
        val, tvals, idxs, deriv, p, continuity::Symbol = :left
    )
    projected = interp(tvals, idxs, deriv, p, continuity)
    val isa RecursiveArrayTools.AbstractDiffEqArray ? copyto!(val.u, projected) :
        copyto!(val, projected)
    return val
end

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::AbstractGLEEAlgorithm, args...;
        kwargs...
    )
    integrator = SciMLBase.__init(prob, alg, args...; kwargs...)
    SciMLBase.solve!(integrator)
    return _split_global_error_solution(integrator.sol, prob)
end

# Turn the extended (y, ε) ArrayPartition solution into an ordinary solution
# over the user's problem: `u` holds the solution partition, `global_error`
# holds the error partition, and the interpolation projects onto the solution.
function _split_global_error_solution(ext_sol, prob)
    us = [u.x[1] for u in ext_sol.u]
    ges = [u.x[2] for u in ext_sol.u]
    sol = @set ext_sol.interp = GLEESolutionInterpolation(ext_sol.interp)
    sol = @set sol.prob = prob
    sol = @set sol.global_error = ges
    # Set `u` last: `setproperties` recomputes the solution's `T`/`N` type
    # parameters from the new (unpartitioned) `u`.
    return @set sol.u = us
end

"""
    global_error_estimate(sol)
    global_error_estimate(sol, i)

Extract the global error estimate from a solution computed with a GLEE method
([`GLEE23`](@ref), [`GLEE24`](@ref), [`GLEE35`](@ref), [`MM5GEE`](@ref)).

These solvers report the global error through the standard SciMLBase interface:
`sol.global_error` is a vector matching `sol.u`, with `sol.global_error[i]` the
estimated global error of `sol.u[i]` at `sol.t[i]` (see
[`SciMLBase.has_global_error`](@ref)). `global_error_estimate(sol)` returns that
vector and `global_error_estimate(sol, i)` its `i`-th element.
"""
function global_error_estimate(sol::SciMLBase.AbstractODESolution)
    sol.global_error === nothing && throw(
        ArgumentError(
            "global_error_estimate expects a solution produced by a GLEE method"
        )
    )
    return sol.global_error
end

function global_error_estimate(sol::SciMLBase.AbstractODESolution, i::Integer)
    return global_error_estimate(sol)[i]
end
