# TODO: Optimize cache size
mutable struct AN5ConstantCache{zType, lType, dtsType, dType, tsit5Type} <:
               OrdinaryDiffEqConstantCache
    # `z` is the Nordsieck vector
    z::zType
    # `l` is used for the corrector iteration
    l::Vector{lType}
    # `m` is a tmp vector that is used for calculating `l`
    m::Vector{lType}
    # `c_LTE` is used for the error estimation for the current order
    c_LTE::lType
    c_conv::lType
    # `dts` stores `dt`s
    dts::dtsType
    # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
    Δ::dType
    # `Tsit5` for the first step
    tsit5tab::tsit5Type
    order::Int
end

function alg_cache(alg::AN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    N = 5
    z = [zero(rate_prototype) for i in 1:(N + 1)]
    Δ = u
    l = fill(zero(tTypeNoUnits), N + 1)
    m = zero(l)
    c_LTE = c_conv = zero(tTypeNoUnits)
    dts = fill(zero(dt), 6)
    tsit5tab = Tsit5ConstantCache()
    AN5ConstantCache(z, l, m, c_LTE, c_conv, dts, Δ, tsit5tab, 1)
end

mutable struct AN5Cache{uType, dType, rateType, zType, lType, dtsType, tsit5Type} <:
               OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    Δ::dType
    # Error estimation
    atmp::dType
    fsalfirst::rateType
    ratetmp::rateType
    # `z` is the Nordsieck vector
    z::zType
    # `l` is used for the corrector iteration
    l::Vector{lType}
    # `m` is a tmp vector that is used for calculating `l`
    m::Vector{lType}
    # `c_LTE` is used for the error estimation for the current order
    c_LTE::lType
    c_conv::lType
    # `dts` stores `dt`s
    dts::dtsType
    # `Tsit5` for the first step
    tsit5cache::tsit5Type
    order::Int
end

function alg_cache(alg::AN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    #################################################
    # Tsit5
    # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    tsit5cache = Tsit5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp,
        trivial_limiter!, trivial_limiter!, False())
    #################################################
    N = 5
    Δ = similar(atmp)
    recursivefill!(Δ, false)
    l = fill(zero(tTypeNoUnits), N + 1)
    m = zero(l)
    c_LTE = c_conv = zero(tTypeNoUnits)
    dts = fill(zero(dt), 6)
    fsalfirst = zero(rate_prototype)
    z = [zero(rate_prototype) for i in 1:(N + 1)]
    for i in 1:(N + 1)
        z[i] = zero(rate_prototype)
    end
    ratetmp = zero(rate_prototype)

    AN5Cache(u, uprev, tmp, Δ, atmp, fsalfirst, ratetmp,
        z, l, m, c_LTE, c_conv, dts,
        tsit5cache, 1)
end

mutable struct JVODEConstantCache{zType, lType, dtsType, dType, tsit5Type, etaType} <:
               OrdinaryDiffEqConstantCache
    # `z` is the Nordsieck vector
    z::zType
    # `l` is used for the corrector iteration
    l::Vector{lType}
    # `m` is a tmp vector that is used for calculating `l`
    m::Vector{lType}
    # `c_LTE₊₁` is used for the error estimation for the current order + 1
    c_LTE₊₁::lType
    # `c_LTE` is used for the error estimation for the current order
    c_LTE::lType
    # `c_LTE₋₁` is used for the error estimation for the current order - 1
    c_LTE₋₁::lType
    # `c_conv` is used in convergence test
    c_conv::lType
    # `c_𝒟` is used to get the order q+2 derivative vector
    c_𝒟::lType
    prev_𝒟::lType
    # `dts` stores `dt`s
    dts::dtsType
    # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
    Δ::dType
    # `Tsit5` for the first step
    tsit5tab::tsit5Type
    L::Int
    # same with `order` or `q`
    order::Int
    nextorder::Int
    # number of steps to take before considering to change order
    n_wait::Int
    # `η` is `dtₙ₊₁/dtₙ`
    η::etaType
    ηold::etaType
    ηq::etaType
    η₊₁::etaType
    η₋₁::etaType
    maxη::etaType
end

function alg_cache(alg::JVODE, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    N = 12
    z = [rate_prototype for i in 1:(N + 1)]
    Δ = u
    l = fill(zero(tTypeNoUnits), N + 1)
    m = zero(l)
    c_LTE₊₁ = c_LTE = c_LTE₋₁ = c_conv = c_𝒟 = prev_𝒟 = zero(tTypeNoUnits)
    dts = fill(zero(dt), N + 1)
    tsit5tab = Tsit5ConstantCache()
    η = zero(dt / dt)
    JVODEConstantCache(z, l, m,
        c_LTE₊₁, c_LTE, c_LTE₋₁, c_conv, c_𝒟, prev_𝒟,
        dts, Δ, tsit5tab, 2, 1, 1, 2, η, η, η, η, η, η)
end

mutable struct JVODECache{
    uType,
    rateType,
    zType,
    lType,
    dtsType,
    dType,
    etaType,
    tsit5Type
} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    fsalfirst::rateType
    ratetmp::rateType
    # `z` is the Nordsieck vector
    z::zType
    # `l` is used for the corrector iteration
    l::Vector{lType}
    # `m` is a tmp vector that is used for calculating `l`
    m::Vector{lType}
    # `c_LTE₊₁` is used for the error estimation for the current order + 1
    c_LTE₊₁::lType
    # `c_LTE` is used for the error estimation for the current order
    c_LTE::lType
    # `c_LTE₋₁` is used for the error estimation for the current order - 1
    c_LTE₋₁::lType
    # `c_conv` is used in convergence test
    c_conv::lType
    # `c_𝒟` is used to get the order q+2 derivative vector
    c_𝒟::lType
    prev_𝒟::lType
    # `dts` stores `dt`s
    dts::dtsType
    # `Δ` is the difference between the predictor `uₙ₀` and `uₙ`
    Δ::dType
    # Error estimation
    atmp::dType
    # `Tsit5` for the first step
    tsit5cache::tsit5Type
    L::Int
    # same with `order` or `q`
    order::Int
    nextorder::Int
    # number of steps to take before considering to change order
    n_wait::Int
    # `η` is `dtₙ₊₁/dtₙ`
    η::etaType
    ηold::etaType
    ηq::etaType
    η₊₁::etaType
    η₋₁::etaType
    maxη::etaType
end

function alg_cache(alg::JVODE, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    #################################################
    # Tsit5
    # Cannot alias pointers, since we have to use `k`s to start the Nordsieck vector
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k6 = zero(rate_prototype)
    k7 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    tsit5cache = Tsit5Cache(u, uprev, k1, k2, k3, k4, k5, k6, k7, utilde, tmp, atmp,
        trivial_limiter!, trivial_limiter!, False())
    #################################################
    fsalfirst = zero(rate_prototype)
    N = 12
    z = [zero(rate_prototype) for i in 1:(N + 1)]
    Δ = similar(u, uEltypeNoUnits)
    recursivefill!(Δ, false)
    l = fill(zero(tTypeNoUnits), N + 1)
    m = zero(l)
    c_LTE₊₁ = c_LTE = c_LTE₋₁ = c_conv = c_𝒟 = prev_𝒟 = zero(tTypeNoUnits)
    dts = fill(zero(dt), N + 1)
    η = zero(dt / dt)
    #################################################
    # Nordsieck Vector
    z[1] = zero(rate_prototype)
    z[2] = zero(rate_prototype)
    z[3] = zero(rate_prototype)
    z[4] = zero(rate_prototype)
    z[5] = zero(rate_prototype)
    z[6] = zero(rate_prototype)
    z[7] = zero(rate_prototype)
    z[8] = zero(rate_prototype)
    z[9] = zero(rate_prototype)
    z[10] = zero(rate_prototype)
    z[11] = zero(rate_prototype)
    z[12] = zero(rate_prototype)
    z[13] = zero(rate_prototype)
    ratetmp = zero(rate_prototype)
    #################################################
    JVODECache(u, uprev, tmp, fsalfirst, ratetmp,
        z, l, m,
        c_LTE₊₁, c_LTE, c_LTE₋₁, c_conv, c_𝒟, prev_𝒟,
        dts, Δ, atmp, tsit5cache, 2, 1, 1, 2, η, η, η, η, η, η)
end

function get_fsalfirstlast(cache::Union{JVODECache, AN5Cache}, u)
    get_fsalfirstlast(cache.tsit5cache, u)
end
