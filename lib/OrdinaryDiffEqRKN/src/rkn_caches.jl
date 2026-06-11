abstract type NystromMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::NystromMutableCache, u) = (cache.fsalfirst, cache.k)

## Generic velocity-independent Nyström caches

struct NystromVIConstantCache{TabType} <: NystromConstantCache
    tab::TabType
end

@cache struct NystromVICache{uType, rateType, reducedRateType, uNoUnitsType, TabType} <:
    NystromMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    ks::Vector{reducedRateType}   # stage derivatives k2..kN (length nstages-1)
    k::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::TabType
end

# DPRKN6 runs through the generic VI caches but is identified by its tableau type so
# its specialized dense-output interpolant can dispatch.
const DPRKN6Caches = Union{
    NystromVIConstantCache{<:DPRKN6Tableau},
    NystromVICache{<:Any, <:Any, <:Any, <:Any, <:DPRKN6Tableau},
}

## Generic velocity-dependent Nyström caches

struct NystromVDConstantCache{T, T2} <: NystromConstantCache
    tab::NystromVDTableau{T, T2}
end

@cache struct NystromVDCache{uType, rateType, reducedRateType, uNoUnitsType, T, T2} <:
    NystromMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    ks::Vector{reducedRateType}   # stage derivatives k2..kN (length nstages-1)
    k::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    tab::NystromVDTableau{T, T2}
end

function alg_cache(
        alg::Nystrom4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = Nystrom4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVDCache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::Nystrom4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVDConstantCache(
        Nystrom4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::FineRKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = FineRKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVDCache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::FineRKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVDConstantCache(
        FineRKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::FineRKN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = FineRKN5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVDCache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::FineRKN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVDConstantCache(
        FineRKN5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

@cache struct Nystrom4VelocityIndependentCache{uType, rateType, reducedRateType} <:
    NystromMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k₂::reducedRateType
    k₃::reducedRateType
    k::rateType
    tmp::uType
end

function alg_cache(
        alg::Nystrom4VelocityIndependent, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = Nystrom4VelocityIndependentTableau(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

struct Nystrom4VelocityIndependentConstantCache <: NystromConstantCache end

function alg_cache(
        alg::Nystrom4VelocityIndependent, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Nystrom4VelocityIndependentTableau(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    return NystromVIConstantCache(tab)
end

@cache struct IRKN3Cache{uType, rateType, TabType} <: NystromMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k₂::rateType
    k::rateType
    tmp::uType
    tmp2::rateType
    onestep_cache::Nystrom4VelocityIndependentCache
    tab::TabType
end

function alg_cache(
        alg::IRKN3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    k₃ = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    tab = IRKN3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return IRKN3Cache(
        u, uprev, uprev2, k₁, k₂, k, tmp, k₃,
        Nystrom4VelocityIndependentCache(u, uprev, k₁, k₂.x[2], k₃.x[2], k, tmp),
        tab
    )
end

function alg_cache(
        alg::IRKN3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return IRKN3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

@cache struct IRKN4Cache{uType, rateType, TabType} <: NystromMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k₂::rateType
    k₃::rateType
    k::rateType
    tmp::uType
    tmp2::rateType
    onestep_cache::Nystrom4VelocityIndependentCache
    tab::TabType
end

function alg_cache(
        alg::IRKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    k₃ = zero(rate_prototype)
    k = zero(rate_prototype)
    tmp = zero(u)
    tmp2 = zero(rate_prototype)
    tab = IRKN4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return IRKN4Cache(
        u, uprev, uprev2, k₁, k₂, k₃, k, tmp, tmp2,
        Nystrom4VelocityIndependentCache(u, uprev, k₁, k₂.x[2], k₃.x[2], k, tmp),
        tab
    )
end

function alg_cache(
        alg::IRKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return IRKN4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function alg_cache(
        alg::Nystrom5VelocityIndependent, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = Nystrom5VelocityIndependentTableau(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::Nystrom5VelocityIndependent, u, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = Nystrom5VelocityIndependentTableau(
        constvalue(uBottomEltypeNoUnits),
        constvalue(tTypeNoUnits)
    )
    return NystromVIConstantCache(tab)
end

function alg_cache(
        alg::DPRKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = DPRKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::DPRKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        DPRKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::DPRKN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = DPRKN5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::DPRKN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        DPRKN5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::DPRKN6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = DPRKN6Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::DPRKN6, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        DPRKN6Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::DPRKN6FM, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = DPRKN6FMTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::DPRKN6FM, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        DPRKN6FMTableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::DPRKN8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = DPRKN8Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::DPRKN8, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        DPRKN8Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::DPRKN12, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = DPRKN12Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::DPRKN12, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        DPRKN12Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::ERKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = ERKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::ERKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        ERKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::ERKN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = ERKN5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::ERKN5, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        ERKN5Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::ERKN7, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = ERKN7Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVICache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::ERKN7, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVIConstantCache(
        ERKN7Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end

function alg_cache(
        alg::RKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    reduced_rate_prototype = rate_prototype.x[2]
    tab = RKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    nstages = length(tab.b)
    k1 = zero(rate_prototype)
    ks = [zero(reduced_rate_prototype) for _ in 2:nstages]
    k = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return NystromVDCache(u, uprev, k1, ks, k, utilde, tmp, atmp, tab)
end

function alg_cache(
        alg::RKN4, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return NystromVDConstantCache(
        RKN4Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    )
end
