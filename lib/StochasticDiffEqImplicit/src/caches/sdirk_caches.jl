@cache mutable struct ImplicitEMCache{uType, rateType, N, noiseRateType, dWType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    gtmp::noiseRateType
    gtmp2::rateType
    dW_cache::dWType
    nlsolver::N
end

function alg_cache(
        alg::ImplicitEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    gtmp = zero(noise_rate_prototype)
    if is_diagonal_noise(prob)
        gtmp2 = gtmp
        dW_cache = nothing
    else
        gtmp2 = zero(rate_prototype)
        dW_cache = zero(ΔW)
    end

    γ, c = alg.theta, zero(t)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)
    return ImplicitEMCache(u, uprev, fsalfirst, gtmp, gtmp2, dW_cache, nlsolver)
end

mutable struct ImplicitEMConstantCache{N} <: StochasticDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::ImplicitEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = alg.theta, zero(t)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return ImplicitEMConstantCache(nlsolver)
end

@cache mutable struct ImplicitEulerHeunCache{uType, rateType, N, noiseRateType, dWType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    gtmp::noiseRateType
    gtmp2::rateType
    gtmp3::noiseRateType
    nlsolver::N
    dW_cache::dWType
end

u_cache(c::ImplicitEulerHeunCache) = (c.uprev2, c.nlsolver.z, c.nlsolver.dz)
du_cache(c::ImplicitEulerHeunCache) = (c.nlsolver.k, c.fsalfirst)

function alg_cache(
        alg::ImplicitEulerHeun, prob, u, ΔW, ΔZ, p,
        rate_prototype, noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    gtmp = zero(noise_rate_prototype)
    gtmp2 = zero(rate_prototype)

    if is_diagonal_noise(prob)
        gtmp3 = gtmp2
        dW_cache = nothing
    else
        gtmp3 = zero(noise_rate_prototype)
        dW_cache = zero(ΔW)
    end

    γ, c = alg.theta, zero(t)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)
    return ImplicitEulerHeunCache(u, uprev, fsalfirst, gtmp, gtmp2, gtmp3, nlsolver, dW_cache)
end

mutable struct ImplicitEulerHeunConstantCache{N} <: StochasticDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::ImplicitEulerHeun, prob, u, ΔW, ΔZ, p,
        rate_prototype, noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = alg.theta, zero(t)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return ImplicitEulerHeunConstantCache(nlsolver)
end

@cache mutable struct ImplicitRKMilCache{uType, rateType, N, noiseRateType} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    gtmp::noiseRateType
    gtmp2::noiseRateType
    gtmp3::noiseRateType
    nlsolver::N
end

function alg_cache(
        alg::ImplicitRKMil, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    gtmp = zero(noise_rate_prototype)
    gtmp2 = zero(rate_prototype)
    gtmp3 = zero(rate_prototype)

    γ, c = alg.theta, zero(t)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)
    return ImplicitRKMilCache(u, uprev, fsalfirst, gtmp, gtmp2, gtmp3, nlsolver)
end

mutable struct ImplicitRKMilConstantCache{N} <: StochasticDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::ImplicitRKMil, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = alg.theta, zero(t)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return ImplicitRKMilConstantCache(nlsolver)
end
