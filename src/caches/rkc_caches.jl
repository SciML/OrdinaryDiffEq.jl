mutable struct ROCK2ConstantCache{T, T2, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{46, Int}
    fp1::SVector{46, T}
    fp2::SVector{46, T}
    recf::Vector{T2}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
    min_stage::Int
    max_stage::Int
end
@cache struct ROCK2Cache{uType, rateType, uNoUnitsType, C <: ROCK2ConstantCache} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(alg::ROCK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ROCK2ConstantCache(constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits),
        u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    ROCK2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::ROCK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ROCK2ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits), u)
end

mutable struct ROCK4ConstantCache{T, T2, T3, T4, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{50, Int}
    fpa::Vector{T}
    fpb::Vector{T2}
    fpbe::Vector{T3}
    recf::Vector{T4}
    zprev::zType
    mdeg::Int
    deg_index::Int
    start::Int
    min_stage::Int
    max_stage::Int
end

@cache struct ROCK4Cache{uType, rateType, uNoUnitsType, C <: ROCK4ConstantCache} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    uᵢ₋₃::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(alg::ROCK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ROCK4ConstantCache(constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits),
        u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    uᵢ₋₃ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    ROCK4Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, uᵢ₋₃, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::ROCK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ROCK4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits), u)
end

mutable struct RKCConstantCache{zType} <: OrdinaryDiffEqConstantCache
    #to match the types to call maxeig!
    zprev::zType
end
@cache struct RKCCache{uType, rateType, uNoUnitsType, C <: RKCConstantCache} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(alg::RKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = RKCConstantCache(u)
    gprev = zero(u)
    gprev2 = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    RKCCache(u, uprev, gprev, gprev2, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::RKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    RKCConstantCache(u)
end

@cache mutable struct IRKCConstantCache{uType, rateType, N} <: OrdinaryDiffEqConstantCache
    minm::Int
    zprev::uType
    nlsolver::N
    du₁::rateType
    du₂::rateType
end

@cache mutable struct IRKCCache{uType, rateType, uNoUnitsType, N, C <: IRKCConstantCache} <:
                      OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    gprev::uType
    gprev2::uType
    fsalfirst::rateType
    f1ⱼ₋₁::rateType
    f1ⱼ₋₂::rateType
    f2ⱼ₋₁::rateType
    atmp::uNoUnitsType
    nlsolver::N
    du₁::rateType
    du₂::rateType
    constantcache::C
end

function alg_cache(alg::IRKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1.0, 1.0
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false))
    zprev = u
    du₁ = rate_prototype
    du₂ = rate_prototype
    IRKCConstantCache(50, zprev, nlsolver, du₁, du₂)
end

function alg_cache(alg::IRKC, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1.0, 1.0
    nlsolver = build_nlsolver(alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true))

    gprev = zero(u)
    gprev2 = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    zprev = zero(u)
    f1ⱼ₋₁ = zero(rate_prototype)
    f1ⱼ₋₂ = zero(rate_prototype)
    f2ⱼ₋₁ = zero(rate_prototype)
    du₁ = zero(rate_prototype)
    du₂ = zero(rate_prototype)
    constantcache = IRKCConstantCache(50, zprev, nlsolver, du₁, du₂)
    IRKCCache(u, uprev, gprev, gprev2, fsalfirst, f1ⱼ₋₁, f1ⱼ₋₂, f2ⱼ₋₁, atmp, nlsolver, du₁,
        du₂, constantcache)
end

mutable struct ESERK4ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{46, Int}
    Cᵤ::SVector{4, Int}
    Cₑ::SVector{4, Int}
    zprev::zType
    Bᵢ::Vector{T}
    mdeg::Int
    start::Int
    internal_deg::Int
end

@cache struct ESERK4Cache{uType, rateType, uNoUnitsType, C <: ESERK4ConstantCache} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(alg::ESERK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ESERK4ConstantCache(u)
    uᵢ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Sᵢ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    ESERK4Cache(u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::ESERK4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ESERK4ConstantCache(u)
end

mutable struct ESERK5ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{49, Int}
    Cᵤ::SVector{5, Int}
    Cₑ::SVector{5, Int}
    zprev::zType
    Bᵢ::Vector{T}
    mdeg::Int
    start::Int
    internal_deg::Int
end

@cache struct ESERK5Cache{uType, rateType, uNoUnitsType, C <: ESERK5ConstantCache} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(alg::ESERK5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = ESERK5ConstantCache(u)
    uᵢ = zero(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Sᵢ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    ESERK5Cache(u, uprev, uᵢ, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::ESERK5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ESERK5ConstantCache(u)
end

mutable struct SERK2ConstantCache{T, zType} <: OrdinaryDiffEqConstantCache
    ms::SVector{11, Int}
    zprev::zType
    Bᵢ::Vector{T}
    mdeg::Int
    start::Int
    internal_deg::Int
end

@cache struct SERK2Cache{uType, rateType, uNoUnitsType, C <: SERK2ConstantCache} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    uᵢ₋₁::uType
    uᵢ₋₂::uType
    Sᵢ::uType
    tmp::uType
    atmp::uNoUnitsType
    fsalfirst::rateType
    k::rateType
    constantcache::C
end

function alg_cache(alg::SERK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    constantcache = SERK2ConstantCache(u)
    uᵢ₋₁ = zero(u)
    uᵢ₋₂ = zero(u)
    Sᵢ = zero(u)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    fsalfirst = zero(rate_prototype)
    k = zero(rate_prototype)
    SERK2Cache(u, uprev, uᵢ₋₁, uᵢ₋₂, Sᵢ, tmp, atmp, fsalfirst, k, constantcache)
end

function alg_cache(alg::SERK2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    SERK2ConstantCache(u)
end
