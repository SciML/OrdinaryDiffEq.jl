mutable struct SROCK1ConstantCache{zType, tType} <: StochasticDiffEqConstantCache
    ms::SVector{10, Int}
    mη::SVector{10, tType}
    zprev::zType
    mdeg::Int
    optimal_η::tType
end
@cache struct SROCK1Cache{uType, rateType, noiseRateType, C <: SROCK1ConstantCache} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    gₘ₋₁::noiseRateType
    gₘ₋₂::noiseRateType
    tmp::uType
    k::rateType
    fsalfirst::rateType
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::SROCK1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SROCK1ConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::SROCK1, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    gₘ₋₁ = zero(noise_rate_prototype)
    gₘ₋₂ = zero(noise_rate_prototype)
    tmp = uᵢ₋₂             # these 3 variables are dummied to use same memory
    fsalfirst = k
    atmp = zero(rate_prototype)
    constantcache = SROCK1ConstantCache{typeof(t)}(u)
    return SROCK1Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, gₘ₋₁, gₘ₋₂, tmp, k, fsalfirst, atmp, constantcache)
end

mutable struct SROCK2ConstantCache{zType, T} <: StochasticDiffEqConstantCache
    ms::SVector{46, Int}
    recf::Vector{T}
    mσ::SVector{46, T}
    mτ::SVector{46, T}
    recf2::Vector{T}
    mα::SVector{46, T}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
end

@cache struct SROCK2Cache{uType, rateType, noiseRateType, T, C <: SROCK2ConstantCache} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uₓ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Gₛ::noiseRateType
    Gₛ₁::noiseRateType
    WikRange::T
    vec_χ::T
    tmp::uType
    k::rateType
    fsalfirst::rateType
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::SROCK2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SROCK2ConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::SROCK2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    uᵢ = zero(u)
    uₓ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Gₛ = zero(noise_rate_prototype)
    Gₛ₁ = zero(noise_rate_prototype)
    if ΔW isa Union{SArray, Number}
        WikRange = copy(ΔW)
        vec_χ = copy(ΔW)
    else
        WikRange = false .* vec(ΔW)
        vec_χ = false .* vec(ΔW)
    end
    tmp = uᵢ₋₂            # these 2 variables are dummied to use same memory
    fsalfirst = k
    atmp = zero(rate_prototype)
    constantcache = SROCK2ConstantCache{typeof(t)}(u)
    return SROCK2Cache(
        u, uprev, uᵢ, uₓ, uᵢ₋₁, uᵢ₋₂, Gₛ, Gₛ₁, WikRange,
        vec_χ, tmp, k, fsalfirst, atmp, constantcache
    )
end

mutable struct SROCKEMConstantCache{zType, tType} <: StochasticDiffEqConstantCache
    ms::SVector{10, Int}
    mη::SVector{10, tType}
    zprev::zType
    mdeg::Int
    optimal_η::tType
end

@cache struct SROCKEMCache{uType, rateType, noiseRateType, T, C <: SROCKEMConstantCache} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Gₛ::noiseRateType
    Gₛ₁::noiseRateType
    tmp::uType
    k::rateType
    fsalfirst::rateType
    WikRange::T
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::SROCKEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SROCKEMConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::SROCKEM, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Gₛ = zero(noise_rate_prototype)
    if (!alg.strong_order_1 || is_diagonal_noise(prob) || ΔW isa Number || length(ΔW) == 1)
        Gₛ₁ = Gₛ
    else
        Gₛ₁ = zero(noise_rate_prototype)
    end
    if ΔW isa Union{SArray, Number}
        WikRange = copy(ΔW)
    else
        WikRange = false .* vec(ΔW)
    end
    tmp = zero(u)             # these 3 variables are dummied to use same memory
    fsalfirst = k
    atmp = zero(rate_prototype)
    constantcache = SROCKEMConstantCache{typeof(t)}(u)
    return SROCKEMCache(
        u, uprev, uᵢ₋₁, uᵢ₋₂, Gₛ, Gₛ₁, tmp, k, fsalfirst, WikRange, atmp, constantcache
    )
end

mutable struct SKSROCKConstantCache{zType, T} <: StochasticDiffEqConstantCache
    mc::Vector{T}
    mα::Vector{T}
    zprev::zType
end

@cache struct SKSROCKCache{
        uType, rateType, noise_rate_prototype, T, C <: SKSROCKConstantCache,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Gₛ::noise_rate_prototype
    tmp::uType
    k::rateType
    fsalfirst::rateType
    WikRange::T
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::SKSROCK, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SKSROCKConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::SKSROCK, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Gₛ = zero(noise_rate_prototype)
    tmp = uᵢ₋₂             # Dummy variables
    fsalfirst = k
    if ΔW isa Union{SArray, Number}
        WikRange = copy(ΔW)
    else
        WikRange = false .* vec(ΔW)
    end
    atmp = zero(rate_prototype)
    constantcache = SKSROCKConstantCache{typeof(t)}(u)
    return SKSROCKCache(u, uprev, uᵢ₋₁, uᵢ₋₂, Gₛ, tmp, k, fsalfirst, WikRange, atmp, constantcache)
end

mutable struct TangXiaoSROCK2ConstantCache{zType, T} <: StochasticDiffEqConstantCache
    ms::SVector{46, Int}
    recf::Vector{T}
    mσ::SVector{46, T}
    mτ::SVector{46, T}
    recf2::Vector{T}
    mα::SVector{5, T}
    mn̂::SVector{5, Int}
    c1::SVector{13, T}
    c2::SVector{13, T}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
    start_mcs::Int
end

@cache struct TangXiaoSROCK2Cache{
        uType, rateType, noiseRateType, C <: TangXiaoSROCK2ConstantCache,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uₓ::uType
    Û₁::uType
    Û₂::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Gₛ::noiseRateType
    Gₛ₁::noiseRateType
    tmp::uType
    k::rateType
    fsalfirst::rateType
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::TangXiaoSROCK2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return TangXiaoSROCK2ConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::TangXiaoSROCK2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    uᵢ = zero(u)
    uₓ = zero(u)
    Û₁ = zero(u)
    Û₂ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Gₛ = zero(noise_rate_prototype)
    if ΔW isa Number || length(ΔW) == 1 || is_diagonal_noise(prob)
        Gₛ₁ = Gₛ
    else
        Gₛ₁ = zero(noise_rate_prototype)
    end
    tmp = uᵢ₋₂            # these 2 variables are dummied to use same memory
    fsalfirst = k
    atmp = zero(rate_prototype)
    constantcache = TangXiaoSROCK2ConstantCache{typeof(t)}(u)
    return TangXiaoSROCK2Cache(
        u, uprev, uᵢ, uₓ, Û₁, Û₂, uᵢ₋₁, uᵢ₋₂, Gₛ,
        Gₛ₁, tmp, k, fsalfirst, atmp, constantcache
    )
end

mutable struct KomBurSROCK2ConstantCache{zType, T} <: StochasticDiffEqConstantCache
    ms::SVector{46, Int}
    recf::Vector{T}
    mσ::SVector{46, T}
    mτ::SVector{46, T}
    mδ::Vector{T}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
end

@cache struct KomBurSROCK2Cache{
        uType, rateType, noiseRateType, T, C <: KomBurSROCK2ConstantCache,
    } <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    utmp::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    k::rateType
    yₛ₋₁::rateType
    yₛ₋₂::rateType
    yₛ₋₃::rateType
    SXₛ₋₁::uType
    SXₛ₋₂::uType
    SXₛ₋₃::uType
    Gₛ::noiseRateType
    Xₛ₋₁::noiseRateType
    Xₛ₋₂::noiseRateType
    Xₛ₋₃::noiseRateType
    vec_χ::T
    tmp::uType
    fsalfirst::rateType
    WikRange::T
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::KomBurSROCK2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return KomBurSROCK2ConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::KomBurSROCK2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utmp = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    k = zero(rate_prototype)
    yₛ₋₁ = k
    yₛ₋₂ = zero(rate_prototype)
    yₛ₋₃ = zero(rate_prototype)
    Xₛ₋₁ = zero(noise_rate_prototype)
    Xₛ₋₂ = zero(noise_rate_prototype)
    Xₛ₋₃ = zero(noise_rate_prototype)
    vec_χ = false .* vec(ΔW)
    WikRange = false .* vec(ΔW)
    if ΔW isa Number || length(ΔW) == 1 || is_diagonal_noise(prob)
        Gₛ = Xₛ₋₁
        SXₛ₋₁ = utmp
        SXₛ₋₂ = utmp
        SXₛ₋₃ = utmp
    else
        Gₛ = zero(noise_rate_prototype)
        SXₛ₋₁ = zero(uprev)
        SXₛ₋₂ = zero(uprev)
        SXₛ₋₃ = zero(uprev)
    end
    tmp = uᵢ₋₂            # these 3 variables are dummied to use same memory
    fsalfirst = k
    atmp = yₛ₋₂
    constantcache = KomBurSROCK2ConstantCache{typeof(t)}(u)
    return KomBurSROCK2Cache(
        u, uprev, utmp, uᵢ₋₁, uᵢ₋₂, k, yₛ₋₁, yₛ₋₂,
        yₛ₋₃, SXₛ₋₁, SXₛ₋₂, SXₛ₋₃, Gₛ, Xₛ₋₁, Xₛ₋₂, Xₛ₋₃, vec_χ,
        tmp, fsalfirst, WikRange, atmp, constantcache
    )
end

mutable struct SROCKC2ConstantCache{zType, T} <: StochasticDiffEqConstantCache
    ms::SVector{46, Int}
    recf::Vector{T}
    mσ::SVector{46, T}
    mτ::SVector{46, T}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
end

@cache struct SROCKC2Cache{uType, rateType, noiseRateType, T, C <: SROCKC2ConstantCache} <:
    StochasticDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Gₛ::noiseRateType
    Gₛ₁::noiseRateType
    WikRange::T
    tmp::uType
    k::rateType
    fsalfirst::rateType
    atmp::rateType
    constantcache::C
end

function alg_cache(
        alg::SROCKC2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SROCKC2ConstantCache{typeof(t)}(u)
end

function alg_cache(
        alg::SROCKC2, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype)
    uᵢ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Gₛ = zero(noise_rate_prototype)
    Gₛ₁ = zero(noise_rate_prototype)
    WikRange = false .* vec(ΔW)
    tmp = zero(u)
    fsalfirst = k           # this variables are dummied to use same memory
    atmp = zero(rate_prototype)
    constantcache = SROCKC2ConstantCache{typeof(t)}(u)
    return SROCKC2Cache(
        u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Gₛ, Gₛ₁, WikRange, tmp, k, fsalfirst, atmp, constantcache
    )
end
