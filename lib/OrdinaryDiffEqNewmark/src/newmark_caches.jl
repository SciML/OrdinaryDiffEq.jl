@cache struct NewmarkBetaCache{uType, rateType, parameterType, N} <:
              OrdinaryDiffEqMutableCache
    u::uType # Current solution
    uprev::uType # Previous solution
    fsalfirst::rateType
    β::parameterType # newmark parameter 1
    γ::parameterType # newmark parameter 2
    nlcache::N # Inner solver
    tmp::uType # temporary, because it is required.
end

function alg_cache(alg::NewmarkBeta, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    (; β, γ) = alg
    fsalfirst = zero(rate_prototype)

    # Temporary terms
    aₙ = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ
    )
    aₙ₊₁ = zero(u.x[1])
    prob = NonlinearProblem{true}(newmark_discretized_residual!, aₙ₊₁, evalcache)
    nlcache = init(prob, alg.nlsolve)

    tmp = zero(u)
    NewmarkBetaCache(u, uprev, fsalfirst, β, γ, nlcache, tmp)
end

@cache struct NewmarkBetaConstantCache{uType, rateType, parameterType, N, N2} <:
              OrdinaryDiffEqConstantCache
    u::uType # Current solution
    uprev::uType # Previous solution
    fsalfirst::rateType
    β::parameterType # newmark parameter 1
    γ::parameterType # newmark parameter 2
    prob::N # Inner solver
    nlsolver::N2
    tmp::uType # temporary, because it is required.
end

function alg_cache(alg::NewmarkBeta, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    (; β, γ) = alg
    fsalfirst = zero(rate_prototype)

    # Temporary terms
    aₙ = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ
    )
    aₙ₊₁ = zero(u.x[1])

    tmp = zero(u)
    NewmarkBetaConstantCache(u, uprev, fsalfirst, β, γ, nothing, alg.nlsolve, tmp)
end

function get_fsalfirstlast(cache::Union{NewmarkBetaCache, NewmarkBetaConstantCache}, u)
    (cache.fsalfirst, u)
end
