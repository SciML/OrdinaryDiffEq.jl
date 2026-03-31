struct TauLeapingConstantCache <: StochasticDiffEqConstantCache end

@cache struct TauLeapingCache{uType, rateType} <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    newrate::rateType
    EEstcache::rateType
end

function alg_cache(
        alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return TauLeapingConstantCache()
end

function alg_cache(
        alg::TauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    newrate = zero(jump_rate_prototype)
    EEstcache = zero(jump_rate_prototype)
    return TauLeapingCache(u, uprev, tmp, newrate, EEstcache)
end

function alg_cache(
        alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return TauLeapingConstantCache()
end

function alg_cache(
        alg::CaoTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(u)
    return TauLeapingCache(u, uprev, tmp, nothing, nothing)
end

# ImplicitTauLeaping cache
# First-order implicit (backward Euler) tau-leaping method
# Uses standard nlsolver infrastructure from OrdinaryDiffEqNonlinearSolve

struct ImplicitTauLeapingConstantCache{rateType, N} <: StochasticDiffEqConstantCache
    poisson_counts::rateType  # Storage for Poisson random variates
    rate_at_uprev::rateType   # a(X_n)
    nlsolver::N
end

@cache struct ImplicitTauLeapingCache{uType, rateType, N} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    poisson_counts::rateType  # k ~ Poisson(dt * a(X_n))
    rate_at_uprev::rateType   # a(X_n)
    nlsolver::N
end

function alg_cache(
        alg::ImplicitTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # γ = 1 (fully implicit), c = 1 (evaluate at t + dt)
    γ, c = one(t), oneunit(t)
    nlsolver = OrdinaryDiffEqNonlinearSolve.build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    poisson_counts = zero(jump_rate_prototype)
    rate_at_uprev = zero(jump_rate_prototype)
    return ImplicitTauLeapingConstantCache(poisson_counts, rate_at_uprev, nlsolver)
end

function alg_cache(
        alg::ImplicitTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # γ = 1 (fully implicit), c = 1 (evaluate at t + dt)
    γ, c = one(t), oneunit(t)
    nlsolver = OrdinaryDiffEqNonlinearSolve.build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    poisson_counts = zero(jump_rate_prototype)
    rate_at_uprev = zero(jump_rate_prototype)
    return ImplicitTauLeapingCache(u, uprev, poisson_counts, rate_at_uprev, nlsolver)
end

# ThetaTrapezoidalTauLeaping cache
# Uses standard nlsolver infrastructure from OrdinaryDiffEqNonlinearSolve
# The nlsolve_f override in integrator_utils.jl provides the tau-leaping drift function

struct ThetaTrapezoidalTauLeapingConstantCache{rateType, T, N} <: StochasticDiffEqConstantCache
    poisson_counts::rateType  # Storage for Poisson random variates
    rate_at_uprev::rateType   # a(X_n)
    theta::T
    nlsolver::N
end

@cache struct ThetaTrapezoidalTauLeapingCache{uType, rateType, T, N} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    poisson_counts::rateType  # k ~ Poisson(dt * a(X_n))
    rate_at_uprev::rateType   # a(X_n)
    theta::T
    nlsolver::N
end

function alg_cache(
        alg::ThetaTrapezoidalTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # γ = theta (implicit weight), c = 1 (evaluate at t + dt)
    γ, c = alg.theta, oneunit(t)
    nlsolver = OrdinaryDiffEqNonlinearSolve.build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    poisson_counts = zero(jump_rate_prototype)
    rate_at_uprev = zero(jump_rate_prototype)
    return ThetaTrapezoidalTauLeapingConstantCache(poisson_counts, rate_at_uprev, alg.theta, nlsolver)
end

function alg_cache(
        alg::ThetaTrapezoidalTauLeaping, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # γ = theta (implicit weight), c = 1 (evaluate at t + dt)
    γ, c = alg.theta, oneunit(t)
    nlsolver = OrdinaryDiffEqNonlinearSolve.build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    poisson_counts = zero(jump_rate_prototype)
    rate_at_uprev = zero(jump_rate_prototype)
    return ThetaTrapezoidalTauLeapingCache(u, uprev, poisson_counts, rate_at_uprev, alg.theta, nlsolver)
end
