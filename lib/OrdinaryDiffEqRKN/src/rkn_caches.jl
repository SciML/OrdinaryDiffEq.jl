abstract type NystromMutableCache <: OrdinaryDiffEqMutableCache end
get_fsalfirstlast(cache::NystromMutableCache, u) = (cache.fsalfirst, cache.k)

## Generic velocity-independent Nyström caches

struct NystromVIConstantCache{TabType} <: NystromConstantCache
    tab::TabType
end

@cache struct NystromVICache{uType, rateType, reducedRateType, TabType, TmpC <: TmpCache} <:
    NystromMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    ks::Vector{reducedRateType}   # stage derivatives k2..kN (length nstages-1)
    k::rateType
    # Unified scratch: `tmp` (stage state), `tmp2` (embedded solution, was
    # `utilde`) and `atmp` (error-norm scaling) — migrated fields, so the array
    # count matches the historical cache exactly. Rate slots are `nothing`.
    tmp_cache::TmpC
    tab::TabType
end

# DPRKN6 runs through the generic VI caches but is identified by its tableau type so
# its specialized dense-output interpolant can dispatch.
const DPRKN6Caches = Union{
    NystromVIConstantCache{<:DPRKN6Tableau},
    NystromVICache{<:Any, <:Any, <:Any, <:DPRKN6Tableau},
}

## Generic velocity-dependent Nyström caches

struct NystromVDConstantCache{T, T2} <: NystromConstantCache
    tab::NystromVDTableau{T, T2}
end

@cache struct NystromVDCache{uType, rateType, reducedRateType, T, T2, TmpC <: TmpCache} <:
    NystromMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    ks::Vector{reducedRateType}   # stage derivatives k2..kN (length nstages-1)
    k::rateType
    # Unified scratch: `tmp` (stage state), `tmp2` (embedded solution, was
    # `utilde`) and `atmp` (error-norm scaling) — migrated fields, so the array
    # count matches the historical cache exactly. Rate slots are `nothing`.
    tmp_cache::TmpC
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). No safe rate donors (stages persist via `ks`),
    # so the rate slots stay `nothing`.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVDCache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). No safe rate donors (stages persist via `ks`),
    # so the rate slots stay `nothing`.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVDCache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). No safe rate donors (stages persist via `ks`),
    # so the rate slots stay `nothing`.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVDCache(u, uprev, k1, ks, k, tmp_cache, tab)
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

@cache struct Nystrom4VelocityIndependentCache{
        uType, rateType, reducedRateType, TmpC <: TmpCache,
    } <: NystromMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    k₂::reducedRateType
    k₃::reducedRateType
    k::rateType
    # Unified scratch: only `tmp` is populated (the former inline `tmp` field);
    # this layout never had `utilde`/`atmp`. `tmp2` shares `tmp`'s type
    # parameter, so it donor-aliases the same `tmp` array (no new arrays);
    # `atmp`/`weight` and the rate slots are `nothing`.
    tmp_cache::TmpC
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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

@cache struct IRKN3Cache{uType, rateType, TabType, TmpC <: TmpCache} <:
    NystromMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k₂::rateType
    k::rateType
    # Unified scratch: only `tmp` is populated (the former inline `tmp` field);
    # the TmpCache `tmp2` slot donor-aliases the same `tmp` array (shared type
    # parameter, no new arrays). The cache field `tmp2` below is NOT scratch —
    # it carries stage history (`k1cache`) across steps, so it must not live in
    # (or donate to) the TmpCache slots, which initdt may scribble on mid-solve.
    tmp_cache::TmpC
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
    # Wrap the existing scratch; the onestep bootstrap cache shares the SAME
    # TmpCache so it aliases the same `tmp` array as before (zero new arrays).
    # `tmp`/`tmp2` share the `uType` type parameter, so a tmp-only layout
    # donor-aliases `tmp2` to the same array (`tmp` is dead between steps and
    # never read by dense output). This is safe for initdt: with `atmp === nothing`
    # its state-reuse path (the only place `tmp` and `tmp2` are live together)
    # can never engage.
    tmp_cache = TmpCache(tmp, tmp, nothing, nothing, nothing, nothing)
    tab = IRKN3ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return IRKN3Cache(
        u, uprev, uprev2, k₁, k₂, k, tmp_cache, k₃,
        Nystrom4VelocityIndependentCache(u, uprev, k₁, k₂.x[2], k₃.x[2], k, tmp_cache),
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

@cache struct IRKN4Cache{uType, rateType, TabType, TmpC <: TmpCache} <:
    NystromMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k₂::rateType
    k₃::rateType
    k::rateType
    # Unified scratch: only `tmp` is populated (the former inline `tmp` field);
    # the TmpCache `tmp2` slot donor-aliases the same `tmp` array (shared type
    # parameter, no new arrays). The cache field `tmp2` below is NOT scratch —
    # it carries stage history (`k1cache`) across steps, so it must not live in
    # (or donate to) the TmpCache slots, which initdt may scribble on mid-solve.
    tmp_cache::TmpC
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
    # Wrap the existing scratch; the onestep bootstrap cache shares the SAME
    # TmpCache so it aliases the same `tmp` array as before (zero new arrays).
    # `tmp`/`tmp2` share the `uType` type parameter, so a tmp-only layout
    # donor-aliases `tmp2` to the same array (`tmp` is dead between steps and
    # never read by dense output). This is safe for initdt: with `atmp === nothing`
    # its state-reuse path (the only place `tmp` and `tmp2` are live together)
    # can never engage.
    tmp_cache = TmpCache(tmp, tmp, nothing, nothing, nothing, nothing)
    tab = IRKN4ConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return IRKN4Cache(
        u, uprev, uprev2, k₁, k₂, k₃, k, tmp_cache, tmp2,
        Nystrom4VelocityIndependentCache(u, uprev, k₁, k₂.x[2], k₃.x[2], k, tmp_cache),
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). The stage arrays `ks` feed dense output (DPRKN6's
    # interpolant) so there are no safe rate donors; the rate slots stay
    # `nothing` (RKN algorithms carry no `preallocate_initdt_buffers` knob).
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVICache(u, uprev, k1, ks, k, tmp_cache, tab)
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
    # `tmp`/`tmp2`/`atmp` replace the former inline `tmp`/`utilde`/`atmp`
    # (net-zero array count). No safe rate donors (stages persist via `ks`),
    # so the rate slots stay `nothing`.
    tmp_cache = build_tmp_cache(u, rate_prototype, uEltypeNoUnits)
    return NystromVDCache(u, uprev, k1, ks, k, tmp_cache, tab)
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
