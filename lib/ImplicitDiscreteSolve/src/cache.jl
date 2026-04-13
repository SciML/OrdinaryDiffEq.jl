struct ImplicitDiscreteState{uType, pType, tType}
    u::uType
    p::pType
    t::tType
end

mutable struct IDSolveCache{uType, cType, thetaType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    z::uType
    nlcache::cType
    Θks::thetaType
end

# u === nothing path: no state to evolve (e.g. MTK systems with only callbacks).
# Bypass NonlinearSolve entirely — there is nothing to nonlinear-solve, and
# constructing a NonlinearProblem with empty u would force the nonlinear
# solver to handle a degenerate case. Returning a cache with `nlcache=nothing`
# is dispatched on in `perform_step!` below to make stepping a no-op success.
function alg_cache(
        alg::IDSolve, ::Nothing, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return IDSolveCache(nothing, nothing, nothing, nothing, uBottomEltypeNoUnits[])
end
function alg_cache(
        alg::IDSolve, ::Nothing, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return IDSolveCache(nothing, nothing, nothing, nothing, uBottomEltypeNoUnits[])
end

function alg_cache(
        alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(zero(u), p, t)
    f_nl = (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t)

    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != length(u))
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            u, state
        )
    else
        NonlinearProblem{isinplace(f)}(f_nl, u, state)
    end

    nlcache = init(prob, alg.nlsolve)

    return IDSolveCache(u, uprev, state.u, nlcache, uBottomEltypeNoUnits[])
end

isdiscretecache(cache::IDSolveCache) = true

function alg_cache(
        alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(zero(u), p, t)
    f_nl = (u_next, p) -> f(u_next, p.u, p.p, p.t)

    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != length(u))
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            u, state
        )
    else
        NonlinearProblem{isinplace(f)}(f_nl, u, state)
    end

    nlcache = init(prob, alg.nlsolve)

    # FIXME Use IDSolveConstantCache?
    return IDSolveCache(u, uprev, state.u, nlcache, uBottomEltypeNoUnits[])
end

get_fsalfirstlast(cache::IDSolveCache, rate_prototype) = (nothing, nothing)
