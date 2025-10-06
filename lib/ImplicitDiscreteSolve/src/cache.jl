struct ImplicitDiscreteState{uType, pType, tType}
    u::uType
    p::pType
    t::tType
end

mutable struct IDSolveCache{uType, cType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    z::uType
    nlcache::cType
end

function alg_cache(alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(zero(u), p, t)
    f_nl = (resid, u_next, p) -> f(resid, u_next, p.u, p.p, p.t)

    u_len = length(u)
    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != u_len)
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            u, state)
    else
        NonlinearProblem{isinplace(f)}(f_nl, u, state)
    end

    nlcache = init(prob, alg.nlsolve)

    IDSolveCache(u, uprev, state.u, nlcache)
end

isdiscretecache(cache::IDSolveCache) = true

# struct IDSolveConstantCache <: OrdinaryDiffEqConstantCache
#     prob::Union{Nothing, SciMLBase.AbstractNonlinearProblem}
# end

function alg_cache(alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    state = ImplicitDiscreteState(zero(u), p, t)
    f_nl = (u_next, p) -> f(u_next, p.u, p.p, p.t)

    u_len = length(u)
    nlls = !isnothing(f.resid_prototype) && (length(f.resid_prototype) != u_len)
    prob = if nlls
        NonlinearLeastSquaresProblem{isinplace(f)}(
            NonlinearFunction(f_nl; resid_prototype = f.resid_prototype),
            u, state)
    else
        NonlinearProblem{isinplace(f)}(f_nl, u, state)
    end

    nlcache = init(prob, alg.nlsolve)

    # FIXME Use IDSolveConstantCache
    IDSolveCache(u, uprev, state.u, nlcache)
end

get_fsalfirstlast(cache::IDSolveCache, rate_prototype) = (nothing, nothing)
