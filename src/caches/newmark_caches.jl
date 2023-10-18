@cache struct NewmarkCache{uType, rateType, parameterType, N} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    β::parameterType
    γ::parameterType
    tmp::uType
    nlsolver::N
end

function alg_cache(alg::Newmark, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck, 
    ::Val{true}) {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}

    β = alg.β
    γ = alg.γ
    tmp = zero(u)

    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))
    
    NewmarkCache(u, uprev, β, γ, tmp, nlsolver)
end
