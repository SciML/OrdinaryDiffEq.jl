abstract type ABMMutableCache <: OrdinaryDiffEqMutableCache end
abstract type ABMVariableCoefficientMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::ABMMutableCache, u) = (cache.fsalfirst, cache.k)
function get_fsalfirstlast(cache::ABMVariableCoefficientMutableCache, u)
    return (cache.fsalfirst, cache.k4)
end
@cache mutable struct AB3Cache{uType, rateType, Thread} <: ABMMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    ralk2::rateType
    k::rateType
    tmp::uType
    step::Int
    thread::Thread
end

@cache mutable struct AB3ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
    k2::rateType
    k3::rateType
    step::Int
end

function alg_cache(
        alg::AB3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    ralk2 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    return AB3Cache(u, uprev, fsalfirst, k2, k3, ralk2, k, tmp, 1, alg.thread)
end

function alg_cache(
        alg::AB3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    k3 = rate_prototype
    return AB3ConstantCache(k2, k3, 1)
end

@cache mutable struct ABM32Cache{uType, rateType, Thread} <: ABMMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    ralk2::rateType
    k::rateType
    tmp::uType
    step::Int
    thread::Thread
end

@cache mutable struct ABM32ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
    k2::rateType
    k3::rateType
    step::Int
end

function alg_cache(
        alg::ABM32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    ralk2 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    return ABM32Cache(u, uprev, fsalfirst, k2, k3, ralk2, k, tmp, 1, alg.thread)
end

function alg_cache(
        alg::ABM32, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    k3 = rate_prototype
    return ABM32ConstantCache(k2, k3, 1)
end

@cache mutable struct AB4Cache{uType, rateType, Thread} <: ABMMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    ralk2::rateType
    k::rateType
    tmp::uType
    t2::rateType
    t3::rateType
    t4::rateType
    step::Int
    thread::Thread
end

@cache mutable struct AB4ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
    k2::rateType
    k3::rateType
    k4::rateType
    step::Int
end

function alg_cache(
        alg::AB4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    ralk2 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    t2 = zero(rate_prototype)
    t3 = zero(rate_prototype)
    t4 = zero(rate_prototype)
    return AB4Cache(u, uprev, fsalfirst, k2, k3, k4, ralk2, k, tmp, t2, t3, t4, 1, alg.thread)
end

function alg_cache(
        alg::AB4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    k3 = rate_prototype
    k4 = rate_prototype
    return AB4ConstantCache(k2, k3, k4, 1)
end

@cache mutable struct ABM43Cache{uType, rateType, Thread} <: ABMMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    ralk2::rateType
    k::rateType
    tmp::uType
    t2::rateType
    t3::rateType
    t4::rateType
    t5::rateType
    t6::rateType
    t7::rateType
    step::Int
    thread::Thread
end

@cache mutable struct ABM43ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
    k2::rateType
    k3::rateType
    k4::rateType
    step::Int
end

function alg_cache(
        alg::ABM43, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    ralk2 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    t2 = zero(rate_prototype)
    t3 = zero(rate_prototype)
    t4 = zero(rate_prototype)
    t5 = zero(rate_prototype)
    t6 = zero(rate_prototype)
    t7 = zero(rate_prototype)
    return ABM43Cache(
        u, uprev, fsalfirst, k2, k3, k4, ralk2, k,
        tmp, t2, t3, t4, t5, t6, t7, 1, alg.thread
    )
end

function alg_cache(
        alg::ABM43, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    k3 = rate_prototype
    k4 = rate_prototype
    return ABM43ConstantCache(k2, k3, k4, 1)
end

@cache mutable struct AB5Cache{uType, rateType, Thread} <: ABMMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k::rateType
    tmp::uType
    t2::rateType
    t3::rateType
    t4::rateType
    step::Int
    thread::Thread
end

@cache mutable struct AB5ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    step::Int
end

function alg_cache(
        alg::AB5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    t2 = zero(rate_prototype)
    t3 = zero(rate_prototype)
    t4 = zero(rate_prototype)
    return AB5Cache(u, uprev, fsalfirst, k2, k3, k4, k5, k, tmp, t2, t3, t4, 1, alg.thread)
end

function alg_cache(
        alg::AB5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    k3 = rate_prototype
    k4 = rate_prototype
    k5 = rate_prototype
    return AB5ConstantCache(k2, k3, k4, k5, 1)
end

@cache mutable struct ABM54Cache{uType, rateType, Thread} <: ABMMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    k::rateType
    tmp::uType
    t2::rateType
    t3::rateType
    t4::rateType
    t5::rateType
    t6::rateType
    t7::rateType
    t8::rateType
    step::Int
    thread::Thread
end

@cache mutable struct ABM54ConstantCache{rateType} <: OrdinaryDiffEqConstantCache
    k2::rateType
    k3::rateType
    k4::rateType
    k5::rateType
    step::Int
end

function alg_cache(
        alg::ABM54, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    k4 = zero(rate_prototype)
    k5 = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    t2 = zero(rate_prototype)
    t3 = zero(rate_prototype)
    t4 = zero(rate_prototype)
    t5 = zero(rate_prototype)
    t6 = zero(rate_prototype)
    t7 = zero(rate_prototype)
    t8 = zero(rate_prototype)
    return ABM54Cache(
        u, uprev, fsalfirst, k2, k3, k4, k5, k, tmp,
        t2, t3, t4, t5, t6, t7, t8, 1, alg.thread
    )
end

function alg_cache(
        alg::ABM54, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k2 = rate_prototype
    k3 = rate_prototype
    k4 = rate_prototype
    k5 = rate_prototype
    return ABM54ConstantCache(k2, k3, k4, k5, 1)
end

@cache mutable struct VCAB3ConstantCache{
        TabType, tArrayType, rArrayType, cArrayType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕstar_nm1::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    order::Int
    tab::TabType
    step::Int
end

@cache mutable struct VCAB3Cache{
        uType, rateType, TabType, bs3Type, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    bs3cache::bs3Type
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    atmp::uNoUnitsType
    tmp::uType
    utilde::uType
    tab::TabType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCAB3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 3)
    c = fill(zero(t), 3, 3)
    g = fill(zero(t), 3)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
    for i in 1:3
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 3)
    order = 3
    tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return VCAB3ConstantCache(dts, c, g, ϕ_n, ϕstar_nm1, ϕstar_n, β, order, tab, 1)
end

function alg_cache(
        alg::VCAB3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    bk1 = zero(rate_prototype)
    bk2 = zero(rate_prototype)
    bk3 = zero(rate_prototype)
    bk4 = zero(rate_prototype)
    butilde = zero(u)
    batmp = similar(u, uEltypeNoUnits)
    recursivefill!(batmp, false)
    btmp = zero(u)
    bs3cache = BS3Cache(
        u, uprev, bk1, bk2, bk3, bk4, butilde, btmp, batmp, tab,
        trivial_limiter!, trivial_limiter!, False()
    )
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 3)
    c = fill(zero(t), 3, 3)
    g = fill(zero(t), 3)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
    for i in 1:3
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 3)
    order = 3
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    utilde = zero(u)
    return VCAB3Cache(
        u, uprev, fsalfirst, bs3cache, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕstar_n, β,
        order, atmp, tmp, utilde, tab, 1, alg.thread
    )
end

@cache mutable struct VCAB4ConstantCache{
        rk4constcache, tArrayType, rArrayType, cArrayType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    ϕstar_nm1::rArrayType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    order::Int
    rk4constcache::rk4constcache
    step::Int
end

@cache mutable struct VCAB4Cache{
        uType, rateType, rk4cacheType, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    rk4cache::rk4cacheType
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    atmp::uNoUnitsType
    tmp::uType
    utilde::uType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCAB4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 4)
    c = fill(zero(t), 4, 4)
    g = fill(zero(t), 4)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
    for i in 1:4
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 4)
    order = 4
    rk4constcache = RK4ConstantCache()
    return VCAB4ConstantCache(ϕstar_nm1, dts, c, g, ϕ_n, ϕstar_n, β, order, rk4constcache, 1)
end

function alg_cache(
        alg::VCAB4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    rk1 = zero(rate_prototype)
    rk2 = zero(rate_prototype)
    rk3 = zero(rate_prototype)
    rk4 = zero(rate_prototype)
    rk = zero(rate_prototype)
    rtmp = zero(u)
    ratmp = similar(u, uEltypeNoUnits)
    recursivefill!(ratmp, false)
    rk4cache = RK4Cache(
        u, uprev, rk1, rk2, rk3, rk4, rk, rtmp, ratmp, trivial_limiter!,
        trivial_limiter!, False()
    )
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 4)
    c = fill(zero(t), 4, 4)
    g = fill(zero(t), 4)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
    for i in 1:4
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 4)
    order = 4
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    utilde = zero(u)
    return VCAB4Cache(
        u, uprev, fsalfirst, rk4cache, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕstar_n, β,
        order, atmp, tmp, utilde, 1, alg.thread
    )
end

# VCAB5

@cache mutable struct VCAB5ConstantCache{
        rk4constcache, tArrayType, rArrayType, cArrayType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    ϕstar_nm1::rArrayType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    order::Int
    rk4constcache::rk4constcache
    step::Int
end

@cache mutable struct VCAB5Cache{
        uType, rateType, rk4cacheType, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    rk4cache::rk4cacheType
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    atmp::uNoUnitsType
    tmp::uType
    utilde::uType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCAB5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 5)
    c = fill(zero(t), 5, 5)
    g = fill(zero(t), 5)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
    for i in 1:5
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 5)
    order = 5
    rk4constcache = RK4ConstantCache()
    return VCAB5ConstantCache(ϕstar_nm1, dts, c, g, ϕ_n, ϕstar_n, β, order, rk4constcache, 1)
end

function alg_cache(
        alg::VCAB5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    rk1 = zero(rate_prototype)
    rk2 = zero(rate_prototype)
    rk3 = zero(rate_prototype)
    rk4 = zero(rate_prototype)
    rk = zero(rate_prototype)
    rtmp = zero(u)
    ratmp = similar(u, uEltypeNoUnits)
    recursivefill!(ratmp, false)
    rk4cache = RK4Cache(
        u, uprev, rk1, rk2, rk3, rk4, rk, rtmp, ratmp, trivial_limiter!,
        trivial_limiter!, False()
    )
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 5)
    c = fill(zero(t), 5, 5)
    g = fill(zero(t), 5)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
    for i in 1:5
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 5)
    order = 5
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    utilde = zero(u)
    return VCAB5Cache(
        u, uprev, fsalfirst, rk4cache, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕstar_n, β,
        order, atmp, tmp, utilde, 1, alg.thread
    )
end

# VCABM3

@cache mutable struct VCABM3ConstantCache{
        TabType, tArrayType, rArrayType, cArrayType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕ_np1::rArrayType
    ϕstar_nm1::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    order::Int
    tab::TabType
    step::Int
end

@cache mutable struct VCABM3Cache{
        uType, rateType, TabType, bs3Type, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    bs3cache::bs3Type
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕ_np1::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    atmp::uNoUnitsType
    tmp::uType
    utilde::uType
    tab::TabType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCABM3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 3)
    c = fill(zero(t), 4, 4)
    g = fill(zero(t), 4)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 4)
    for i in 1:3
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 3)
    order = 3
    tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return VCABM3ConstantCache(dts, c, g, ϕ_n, ϕ_np1, ϕstar_nm1, ϕstar_n, β, order, tab, 1)
end

function alg_cache(
        alg::VCABM3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = BS3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    bk1 = zero(rate_prototype)
    bk2 = zero(rate_prototype)
    bk3 = zero(rate_prototype)
    bk4 = zero(rate_prototype)
    butilde = zero(u)
    batmp = similar(u, uEltypeNoUnits)
    recursivefill!(batmp, false)
    btmp = zero(u)
    bs3cache = BS3Cache(
        u, uprev, bk1, bk2, bk3, bk4, butilde, btmp, batmp, tab,
        trivial_limiter!, trivial_limiter!, False()
    )
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 3)
    c = fill(zero(t), 4, 4)
    g = fill(zero(t), 4)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 3)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 3)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 4)
    for i in 1:3
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    for i in 1:4
        ϕ_np1[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 3)
    order = 3
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    utilde = zero(u)
    return VCABM3Cache(
        u, uprev, fsalfirst, bs3cache, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1,
        ϕstar_n, β, order, atmp, tmp, utilde, tab, 1, alg.thread
    )
end

# VCABM4

@cache mutable struct VCABM4ConstantCache{
        rk4constcache, tArrayType, rArrayType, cArrayType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    ϕstar_nm1::rArrayType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕ_np1::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    order::Int
    rk4constcache::rk4constcache
    step::Int
end

@cache mutable struct VCABM4Cache{
        uType, rateType, rk4cacheType, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    rk4cache::rk4cacheType
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕ_np1::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    atmp::uNoUnitsType
    tmp::uType
    utilde::uType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCABM4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 4)
    c = fill(zero(t), 5, 5)
    g = fill(zero(t), 5)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 5)
    for i in 1:4
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 4)
    order = 4
    rk4constcache = RK4ConstantCache()
    return VCABM4ConstantCache(
        ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1, ϕstar_n, β, order, rk4constcache,
        1
    )
end

function alg_cache(
        alg::VCABM4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    rk1 = zero(rate_prototype)
    rk2 = zero(rate_prototype)
    rk3 = zero(rate_prototype)
    rk4 = zero(rate_prototype)
    rk = zero(rate_prototype)
    rtmp = zero(u)
    ratmp = similar(u, uEltypeNoUnits)
    recursivefill!(ratmp, false)
    rk4cache = RK4Cache(
        u, uprev, rk1, rk2, rk3, rk4, rk, rtmp, ratmp, trivial_limiter!,
        trivial_limiter!, False()
    )
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 4)
    c = fill(zero(t), 5, 5)
    g = fill(zero(t), 5)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 4)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 4)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 5)
    for i in 1:4
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    for i in 1:5
        ϕ_np1[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 4)
    order = 4
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    utilde = zero(u)
    return VCABM4Cache(
        u, uprev, fsalfirst, rk4cache, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1,
        ϕstar_n, β, order, atmp, tmp, utilde, 1, alg.thread
    )
end

# VCABM5

@cache mutable struct VCABM5ConstantCache{
        rk4constcache, tArrayType, rArrayType, cArrayType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    ϕstar_nm1::rArrayType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕ_np1::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    order::Int
    rk4constcache::rk4constcache
    step::Int
end

@cache mutable struct VCABM5Cache{
        uType, rateType, rk4cacheType, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    rk4cache::rk4cacheType
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕ_np1::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    atmp::uNoUnitsType
    tmp::uType
    utilde::uType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCABM5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(t), 5)
    c = fill(zero(t), 6, 6)
    g = fill(zero(t), 6)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 6)
    for i in 1:5
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 5)
    order = 5
    rk4constcache = RK4ConstantCache()
    return VCABM5ConstantCache(
        ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1, ϕstar_n, β, order, rk4constcache,
        1
    )
end

function alg_cache(
        alg::VCABM5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    rk1 = zero(rate_prototype)
    rk2 = zero(rate_prototype)
    rk3 = zero(rate_prototype)
    rk4 = zero(rate_prototype)
    rk = zero(rate_prototype)
    rtmp = zero(u)
    ratmp = similar(u, uEltypeNoUnits)
    recursivefill!(ratmp, false)
    rk4cache = RK4Cache(
        u, uprev, rk1, rk2, rk3, rk4, rk, rtmp, ratmp, trivial_limiter!,
        trivial_limiter!, False()
    )
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 5)
    c = fill(zero(t), 6, 6)
    g = fill(zero(t), 6)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 5)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 5)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 6)
    for i in 1:5
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    for i in 1:6
        ϕ_np1[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 5)
    order = 5
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    utilde = zero(u)
    return VCABM5Cache(
        u, uprev, fsalfirst, rk4cache, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1,
        ϕstar_n, β, order, atmp, tmp, utilde, 1, alg.thread
    )
end

# VCABM

@cache mutable struct VCABMConstantCache{
        tArrayType, rArrayType, cArrayType, dtType,
        dtArrayType,
    } <: OrdinaryDiffEqConstantCache
    ϕstar_nm1::rArrayType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::rArrayType
    ϕ_np1::rArrayType
    ϕstar_n::rArrayType
    β::tArrayType
    ξ::dtType
    ξ0::dtType
    order::Int
    max_order::Int
    step::Int
end

@cache mutable struct VCABMCache{
        uType, rateType, dtType, tArrayType, cArrayType,
        uNoUnitsType, coefType, dtArrayType, Thread,
    } <:
    ABMVariableCoefficientMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k4::rateType
    ϕstar_nm1::coefType
    dts::dtArrayType
    c::cArrayType
    g::tArrayType
    ϕ_n::coefType
    ϕ_np1::coefType
    ϕstar_n::coefType
    β::tArrayType
    order::Int
    max_order::Int
    atmp::uNoUnitsType
    tmp::uType
    ξ::dtType
    ξ0::dtType
    utilde::uType
    utildem1::uType
    utildem2::uType
    utildep1::uType
    atmpm1::uNoUnitsType
    atmpm2::uNoUnitsType
    atmpp1::uNoUnitsType
    step::Int
    thread::Thread
end

function alg_cache(
        alg::VCABM, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dts = fill(zero(dt), 13)
    c = fill(zero(t), 13, 13)
    g = fill(zero(t), 13)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 13)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 13)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 13)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 14)
    for i in 1:13
        ϕ_n[i] = copy(rate_prototype)
        ϕstar_nm1[i] = copy(rate_prototype)
        ϕstar_n[i] = copy(rate_prototype)
    end
    β = fill(zero(t), 13)
    ξ = zero(dt)
    ξ0 = zero(dt)
    order = 1
    max_order = 12
    return VCABMConstantCache(
        ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1, ϕstar_n, β, ξ, ξ0, order,
        max_order, 1
    )
end

function alg_cache(
        alg::VCABM, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    fsalfirst = zero(rate_prototype)
    k4 = zero(rate_prototype)
    dts = fill(zero(dt), 13)
    c = fill(zero(t), 13, 13)
    g = fill(zero(t), 13)
    ϕ_n = Vector{typeof(rate_prototype)}(undef, 13)
    ϕstar_nm1 = Vector{typeof(rate_prototype)}(undef, 13)
    ϕstar_n = Vector{typeof(rate_prototype)}(undef, 13)
    ϕ_np1 = Vector{typeof(rate_prototype)}(undef, 14)
    for i in 1:13
        ϕ_n[i] = zero(rate_prototype)
        ϕstar_nm1[i] = zero(rate_prototype)
        ϕstar_n[i] = zero(rate_prototype)
    end
    for i in 1:14
        ϕ_np1[i] = zero(rate_prototype)
    end
    β = fill(zero(t), 13)
    order = 1
    max_order = 12
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    ξ = zero(dt)
    ξ0 = zero(dt)
    utilde = zero(u)
    utildem2 = zero(u)
    utildem1 = zero(u)
    utildep1 = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    atmpm1 = similar(u, uEltypeNoUnits)
    recursivefill!(atmpm1, false)
    atmpm2 = similar(u, uEltypeNoUnits)
    recursivefill!(atmpm2, false)
    atmpp1 = similar(u, uEltypeNoUnits)
    recursivefill!(atmpp1, false)
    return VCABMCache(
        u, uprev, fsalfirst, k4, ϕstar_nm1, dts, c, g, ϕ_n, ϕ_np1, ϕstar_n, β, order,
        max_order, atmp, tmp, ξ, ξ0, utilde, utildem1, utildem2, utildep1, atmpm1,
        atmpm2, atmpp1, 1, alg.thread
    )
end
