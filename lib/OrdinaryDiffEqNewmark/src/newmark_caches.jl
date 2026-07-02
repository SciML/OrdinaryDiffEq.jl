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
        zero(β), zero(β),   # αm = αf = 0 recovers Newmark
        aₙ, vₙ, uₙ,
        DiffCache(copy(aₙ)),
        DiffCache(copy(vₙ)),
        DiffCache(copy(uₙ)),
    )
    aₙ₊₁ = zero(u.x[1])
    prob = NonlinearProblem{true}(discretized_residual!, aₙ₊₁, evalcache)
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
        zero(β), zero(β),   # αm = αf = 0 recovers Newmark
        aₙ, vₙ, uₙ,
        nothing, nothing, nothing,
    )

    tmp = zero(u)
    atmp = zero(u)
    return NewmarkBetaConstantCache(u, uprev, fsalfirst, β, γ, alg.nlsolve, tmp, atmp, thread)
end

@cache struct GeneralizedAlphaCache{uType, rateType, parameterType, N, Thread} <:
    OrdinaryDiffEqMutableCache
    u::uType # Current solution
    uprev::uType # Previous solution
    fsalfirst::rateType
    αm::parameterType # generalized alpha parameter 1
    αf::parameterType # generalized alpha parameter 2
    β::parameterType # newmark parameter 1
    γ::parameterType # newmark parameter 2
    nlcache::N
    tmp::uType
    atmp::uType
    thread::Thread
end

function alg_cache(
        alg::GeneralizedAlpha, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    (; αm, αf, β, γ, thread) = alg
    fsalfirst = zero(rate_prototype)

    aₙ = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        αm, αf,
        aₙ, vₙ, uₙ,
        DiffCache(copy(aₙ)),
        DiffCache(copy(vₙ)),
        DiffCache(copy(uₙ)),
    )
    aₙ₊₁ = zero(u.x[1])
    prob = NonlinearProblem{true}(discretized_residual!, aₙ₊₁, evalcache)
    nlcache = init(prob, alg.nlsolve)

    tmp = zero(u)
    atmp = zero(u)
    return GeneralizedAlphaCache(u, uprev, fsalfirst, αm, αf, β, γ, nlcache, tmp, atmp, thread)
end

@cache struct GeneralizedAlphaConstantCache{uType, rateType, parameterType, N, Thread} <:
    OrdinaryDiffEqConstantCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    αm::parameterType
    αf::parameterType
    β::parameterType
    γ::parameterType
    nlsolver::N
    tmp::uType
    atmp::uType
    thread::Thread
end

function alg_cache(
        alg::GeneralizedAlpha, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    (; αm, αf, β, γ, thread) = alg
    fsalfirst = zero(rate_prototype)

    aₙ = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        αm, αf,
        aₙ, vₙ, uₙ,
        nothing, nothing, nothing,
    )

    tmp = zero(u)
    atmp = zero(u)
    return GeneralizedAlphaConstantCache(u, uprev, fsalfirst, αm, αf, β, γ, alg.nlsolve, tmp, atmp, thread)
end

function get_fsalfirstlast(
        cache::Union{
            NewmarkBetaCache, NewmarkBetaConstantCache,
            GeneralizedAlphaCache, GeneralizedAlphaConstantCache,
        }, u
    )
    return (cache.fsalfirst, cache.fsalfirst)
end
