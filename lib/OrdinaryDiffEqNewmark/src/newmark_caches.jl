@cache struct NewmarkBetaCache{uType, rateType, parameterType, N, Thread} <:
    OrdinaryDiffEqMutableCache
    u::uType # Current solution
    uprev::uType # Previous solution
    fsalfirst::rateType
    β::parameterType # newmark parameter 1
    γ::parameterType # newmark parameter 2
    nlcache::N # Inner solver
    tmp::uType # temporary, because it is required.
    atmp::uType
    thread::Thread
end

function alg_cache(
        alg::NewmarkBeta, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    (; β, γ, thread) = alg
    fsalfirst = zero(rate_prototype)

    # Temporary terms
    aₙ = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ,
        DiffCache(copy(aₙ)),
        DiffCache(copy(vₙ)),
        DiffCache(copy(uₙ)),
    )
    aₙ₊₁ = zero(u.x[1])
    prob = NonlinearProblem{true}(newmark_discretized_residual!, aₙ₊₁, evalcache)
    nlcache = init(prob, alg.nlsolve)

    tmp = zero(u)
    atmp = zero(u)
    return NewmarkBetaCache(u, uprev, fsalfirst, β, γ, nlcache, tmp, atmp, thread)
end

@cache struct NewmarkBetaConstantCache{uType, rateType, parameterType, N, Thread} <:
    OrdinaryDiffEqConstantCache
    u::uType # Current solution
    uprev::uType # Previous solution
    fsalfirst::rateType
    β::parameterType # newmark parameter 1
    γ::parameterType # newmark parameter 2
    nlsolver::N
    tmp::uType # temporary, because it is required.
    atmp::uType
    thread::Thread
end

function alg_cache(
        alg::NewmarkBeta, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    (; β, γ, thread) = alg
    fsalfirst = zero(rate_prototype)

    # Temporary terms
    aₙ = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ,
        nothing, nothing, nothing,
    )

    tmp = zero(u)
    atmp = zero(u)
    return NewmarkBetaConstantCache(u, uprev, fsalfirst, β, γ, alg.nlsolve, tmp, atmp, thread)
end

function get_fsalfirstlast(cache::Union{NewmarkBetaCache, NewmarkBetaConstantCache}, u)
    return (cache.fsalfirst, cache.fsalfirst)
end
