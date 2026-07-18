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
    return SciMLBase.remake(prob; f = extended_f, u0 = u0)
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

function SciMLBase.__solve(
        prob::SciMLBase.AbstractODEProblem, alg::AbstractGLEEAlgorithm, args...;
        kwargs...
    )
    integrator = SciMLBase.__init(prob, alg, args...; kwargs...)
    SciMLBase.solve!(integrator)
    return integrator.sol
end

"""
    global_error_estimate(sol)
    global_error_estimate(sol, i)

Extract the global error estimate from a solution computed with a GLEE method
([`GLEE23`](@ref), [`GLEE24`](@ref), [`GLEE35`](@ref)).

GLEE solutions carry the partitioned state `(y, ε)`: `sol.u[i].x[1]` is the
solution value and `sol.u[i].x[2]` the estimate of its global error at
`sol.t[i]`. `global_error_estimate(sol, i)` returns the error estimate at the
`i`-th saved time point and `global_error_estimate(sol)` returns the vector of
estimates at every saved time point.
"""
function global_error_estimate(sol::SciMLBase.AbstractODESolution)
    return [global_error_estimate(sol, i) for i in eachindex(sol.u)]
end

function global_error_estimate(sol::SciMLBase.AbstractODESolution, i::Integer)
    u = sol.u[i]
    u isa RecursiveArrayTools.ArrayPartition && length(u.x) == 2 || throw(
        ArgumentError(
            "global_error_estimate expects a solution produced by a GLEE method"
        )
    )
    return u.x[2]
end
