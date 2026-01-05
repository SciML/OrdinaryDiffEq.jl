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

function alg_cache(
        alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(isnothing(u) ? nothing : zero(u), p, t)
    f_nl = (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t)

    u_len = isnothing(u) ? 0 : length(u)
    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != u_len)
    unl = isnothing(u) ? Float64[] : u # FIXME nonlinear solve cannot handle nothing for u
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            unl, state
        )
    else
        NonlinearProblem{isinplace(f)}(f_nl, unl, state)
    end

    nlcache = init(prob, alg.nlsolve)

    return IDSolveCache(u, uprev, state.u, nlcache, uBottomEltypeNoUnits[])
end

isdiscretecache(cache::IDSolveCache) = true

function alg_cache(
        alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    @assert !isnothing(u) "Empty u not supported with out of place functions yet."

    state = ImplicitDiscreteState(isnothing(u) ? nothing : zero(u), p, t)
    f_nl = (u_next, p) -> f(u_next, p.u, p.p, p.t)

    u_len = isnothing(u) ? 0 : length(u)
    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != u_len)
    unl = isnothing(u) ? Float64[] : u # FIXME nonlinear solve cannot handle nothing for u
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            unl, state
        )
    else
        NonlinearProblem{isinplace(f)}(f_nl, unl, state)
    end

    nlcache = init(prob, alg.nlsolve)

    # FIXME Use IDSolveConstantCache?
    return IDSolveCache(u, uprev, state.u, nlcache, uBottomEltypeNoUnits[])
end

get_fsalfirstlast(cache::IDSolveCache, rate_prototype) = (nothing, nothing)
