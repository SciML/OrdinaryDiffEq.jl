@cache struct NewmarkBetaCache{uType, rateType, parameterType, N} <: OrdinaryDiffEqMutableCache
    u::uType # Current solution
    uprev::uType # Previous solution
    upred::uType # Predictor solution
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
    upred = zero(u)
    fsalfirst = zero(rate_prototype)

    # Temporary terms
    aₙ     = fsalfirst.x[1]
    vₙ, uₙ = uprev.x
    evalcache = NewmarkDiscretizationCache(
        f, t, p,
        dt, β, γ,
        aₙ, vₙ, uₙ,
    )
    aₙ₊₁   = zero(u.x[1])
    prob = NonlinearProblem{true}(newmark_discretized_residual!, aₙ₊₁, evalcache)
    nlcache = init(prob, alg.nlsolve)

    tmp = zero(u)
    NewmarkBetaCache(u, uprev, upred, fsalfirst, β, γ, nlcache, tmp)
end

function alg_cache(alg::NewmarkBeta, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    @assert false "Using out of place eval with Newmark methods is not supported yet."

    (; β, γ) = alg
    upred = zero(u)
    fsalfirst = zero(rate_prototype)

    @assert 0.0 ≤ β ≤ 0.5
    @assert 0.0 ≤ γ ≤ 1.0

    # ...
    nlsolver = alg.nlsolve

    tmp = zero(u)
    NewmarkBetaCache(u, uprev, upred, fsalfirst, β, γ, nlsolver, tmp)
end

get_fsalfirstlast(cache::NewmarkBetaCache, u) = (cache.fsalfirst, u) # FIXME use fsal
