struct IIF1MConstantCache{vecuType, rhsType, nl_rhsType} <: StochasticDiffEqConstantCache
    uhold::vecuType
    rhs::rhsType
    nl_rhs::nl_rhsType
end

@cache struct IIF1MCache{
        uType, vecuType, DiffCacheType, rhsType, nl_rhsType, rateType,
        rateNoiseType, rateNoiseCollectionType, NoiseTmpType,
    } <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhold::vecuType
    dual_cache::DiffCacheType
    tmp::uType
    rhs::rhsType
    nl_rhs::nl_rhsType
    rtmp1::rateType
    rtmp2::rateNoiseType
    rtmp3::rateNoiseCollectionType
    noise_tmp::NoiseTmpType
end

function alg_cache(
        alg::IIF1M, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uhold = Vector{typeof(u)}(undef, 1)
    tmp = zero(rate_prototype)
    rhs = RHS_IIF1M_Scalar(f, t, t, tmp, p)
    nl_rhs = alg.nlsolve(Val{:init}, rhs, uhold)
    return IIF1MConstantCache(uhold, rhs, nl_rhs)
end

function alg_cache(
        alg::IIF1M, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = similar(u, axes(u))
    rtmp1 = zero(rate_prototype)
    dual_cache = DiffCache(u, Val{determine_chunksize(u, get_chunksize(alg.nlsolve))})
    uhold = vec(u) # this makes uhold the same values as integrator.u
    rhs = RHS_IIF1(f, tmp, t, t, dual_cache, size(u), p)
    nl_rhs = alg.nlsolve(Val{:init}, rhs, uhold)
    noise_tmp = tmp

    rtmp2 = zero(noise_rate_prototype)
    if is_diagonal_noise(prob)
        rtmp3 = rtmp2
    else
        rtmp3 = zero(rate_prototype)
    end
    return IIF1MCache(
        u, uprev, uhold, dual_cache, tmp, rhs, nl_rhs, rtmp1, rtmp2, rtmp3, noise_tmp
    )
end

struct IIF2MConstantCache{vecuType, rhsType, nl_rhsType} <: StochasticDiffEqConstantCache
    uhold::vecuType
    rhs::rhsType
    nl_rhs::nl_rhsType
end

@cache struct IIF2MCache{
        uType, vecuType, DiffCacheType, rhsType, nl_rhsType, rateType,
        rateNoiseType, rateNoiseCollectionType, NoiseTmpType,
    } <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhold::vecuType
    dual_cache::DiffCacheType
    tmp::uType
    rhs::rhsType
    nl_rhs::nl_rhsType
    rtmp1::rateType
    rtmp2::rateNoiseType
    rtmp3::rateNoiseCollectionType
    noise_tmp::NoiseTmpType
end

function alg_cache(
        alg::IIF2M, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uhold = Vector{typeof(u)}(undef, 1)
    tmp = zero(rate_prototype)
    rhs = RHS_IIF2M_Scalar(f, t, t, tmp, p)
    nl_rhs = alg.nlsolve(Val{:init}, rhs, uhold)
    return IIF2MConstantCache(uhold, rhs, nl_rhs)
end

function alg_cache(
        alg::IIF2M, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = similar(u, axes(u))
    rtmp1 = zero(rate_prototype)
    dual_cache = DiffCache(u, Val{determine_chunksize(u, get_chunksize(alg.nlsolve))})
    uhold = vec(u) # this makes uhold the same values as integrator.u
    rhs = RHS_IIF2(f, tmp, t, t, dual_cache, size(u), p)
    nl_rhs = alg.nlsolve(Val{:init}, rhs, uhold)
    noise_tmp = tmp

    rtmp2 = zero(noise_rate_prototype)
    if is_diagonal_noise(prob)
        rtmp3 = rtmp2
    else
        rtmp3 = zero(rate_prototype)
    end
    return IIF2MCache(
        u, uprev, uhold, dual_cache, tmp, rhs, nl_rhs, rtmp1, rtmp2, rtmp3, noise_tmp
    )
end

struct IIF1MilConstantCache{vecuType, rhsType, nl_rhsType} <: StochasticDiffEqConstantCache
    uhold::vecuType
    rhs::rhsType
    nl_rhs::nl_rhsType
end

@cache struct IIF1MilCache{
        uType, vecuType, DiffCacheType, rhsType, nl_rhsType, rateType,
        rateNoiseType, rateNoiseCollectionType, NoiseTmpType,
    } <: StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uhold::vecuType
    dual_cache::DiffCacheType
    tmp::uType
    rhs::rhsType
    nl_rhs::nl_rhsType
    rtmp1::rateType
    rtmp2::rateNoiseType
    rtmp3::rateNoiseCollectionType
    noise_tmp::NoiseTmpType
    gtmp::rateNoiseType
    gtmp2::rateNoiseType
end

function alg_cache(
        alg::IIF1Mil, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    uhold = Vector{typeof(u)}(undef, 1)
    C = zero(rate_prototype)
    rhs = RHS_IIF1M_Scalar(f, t, t, C, p)
    nl_rhs = alg.nlsolve(Val{:init}, rhs, uhold)
    return IIF1MilConstantCache(uhold, rhs, nl_rhs)
end

function alg_cache(
        alg::IIF1Mil, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = similar(u, axes(u))
    rtmp1 = zero(rate_prototype)
    dual_cache = DiffCache(u, Val{determine_chunksize(u, get_chunksize(alg.nlsolve))})
    uhold = vec(u) # this makes uhold the same values as integrator.u
    rhs = RHS_IIF1(f, tmp, t, t, dual_cache, size(u), p)
    nl_rhs = alg.nlsolve(Val{:init}, rhs, uhold)
    noise_tmp = zero(noise_rate_prototype)
    gtmp = zero(noise_rate_prototype)
    gtmp2 = zero(noise_rate_prototype)
    rtmp2 = zero(noise_rate_prototype)
    if is_diagonal_noise(prob)
        rtmp3 = rtmp2
    else
        rtmp3 = zero(rate_prototype)
    end
    return IIF1MilCache(
        u, uprev, uhold, dual_cache, tmp, rhs, nl_rhs,
        rtmp1, rtmp2, rtmp3, noise_tmp, gtmp, gtmp2
    )
end
