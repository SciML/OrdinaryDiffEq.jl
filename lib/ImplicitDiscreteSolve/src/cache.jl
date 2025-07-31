mutable struct ImplicitDiscreteState{uType, pType, tType}
    u::uType
    p::pType
    t_next::tType
end

mutable struct IDSolveCache{uType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    state::ImplicitDiscreteState
    prob::Union{Nothing, SciMLBase.AbstractNonlinearProblem}
end

function alg_cache(alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}

    state = ImplicitDiscreteState(isnothing(u) ? nothing : zero(u), p, t)
    IDSolveCache(u, uprev, state, nothing)
end

isdiscretecache(cache::IDSolveCache) = true

struct IDSolveConstantCache <: OrdinaryDiffEqConstantCache 
    prob::Union{Nothing, SciMLBase.AbstractNonlinearProblem}
end

function alg_cache(alg::IDSolve, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}

    state = ImplicitDiscreteState(isnothing(u) ? nothing : zero(u), p, t)
    IDSolveCache(u, uprev, state, nothing)
end

get_fsalfirstlast(cache::IDSolveCache, rate_prototype) = (nothing, nothing)
