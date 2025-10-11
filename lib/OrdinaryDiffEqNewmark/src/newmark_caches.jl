@cache struct NewmarkBetaCache{uType, rateType, parameterType, N} <: OrdinaryDiffEqMutableCache
    u::uType # Current solution
    uprev::uType # Previous solution
    upred::uType # Predictor solution
    fsalfirst::rateType
    β::parameterType # newmark parameter 1
    γ::parameterType # newmark parameter 2
    nlsolver::N # Inner solver
    tmp::uType # temporary, because it is required.
end

function alg_cache(alg::NewmarkBeta, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}

    β = alg.β
    γ = alg.γ
    upred = zero(u)
    fsalfirst = zero(rate_prototype)

    @assert 0.0 ≤ β ≤ 0.5
    @assert 0.0 ≤ γ ≤ 1.0

    c = 1.0
    γ̂ = NaN # FIXME
    # nlsolver = build_nlsolver(alg, u.x[1], uprev.x[1], p, t, dt, f.f1, rate_prototype.x[1], uEltypeNoUnits,
    #     uBottomEltypeNoUnits, tTypeNoUnits, γ̂, c, Val(true))
    nlsolver = alg.nlsolve

    tmp = zero(u)
    NewmarkBetaCache(u, uprev, upred, fsalfirst, β, γ, nlsolver, tmp)
end

get_fsalfirstlast(cache::NewmarkBetaCache, u) = (cache.fsalfirst, u) # FIXME use fsal
