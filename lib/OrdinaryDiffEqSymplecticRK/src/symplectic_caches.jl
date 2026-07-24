abstract type HamiltonMutableCache <: OrdinaryDiffEqMutableCache end
abstract type HamiltonConstantCache <: OrdinaryDiffEqConstantCache end

@cache struct SymplecticEulerCache{uType, rateType, TmpC <: TmpCache} <: HamiltonMutableCache
    u::uType
    uprev::uType
    # `TmpCache{uType, Nothing, Nothing, Nothing}`: non-adaptive, so no
    # `atmp`/`weight`; no `utilde`, so `tmp2` is opted out. The symplectic
    # algorithm structs are plain singletons (no `preallocate_initdt_buffers`
    # field), so the rate slots are permanently opted out as well — `k` and
    # `fsalfirst` feed the interpolation `k`s and are not legal donors.
    tmp_cache::TmpC
    k::rateType
    fsalfirst::rateType
end

function alg_cache(
        alg::SymplecticEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # The migrated `tmp` (was an inline field) is the only scratch this cache
    # carries; the array count matches the historical cache exactly.
    tmp_cache = TmpCache(zero(u), nothing, nothing, nothing, nothing, nothing)
    return SymplecticEulerCache(u, uprev, tmp_cache, zero(rate_prototype), zero(rate_prototype))
end

struct SymplecticEulerConstantCache <: HamiltonConstantCache end

function alg_cache(
        alg::SymplecticEuler, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return SymplecticEulerConstantCache()
end

@cache struct VelocityVerletCache{uType, rateType, uEltypeNoUnits, TmpC <: TmpCache} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    # `TmpCache{*, Nothing, Nothing, Nothing}` — see SymplecticEulerCache for
    # the slot layout rationale (shared by all caches in this file).
    tmp_cache::TmpC
    k::rateType
    fsalfirst::rateType
    half::uEltypeNoUnits
end

struct VelocityVerletConstantCache{uEltypeNoUnits} <: HamiltonConstantCache
    half::uEltypeNoUnits
end

function alg_cache(
        alg::VelocityVerlet, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(rate_prototype)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    half = uEltypeNoUnits(1 // 2)
    # Migrated `tmp` only (was inline); array count unchanged.
    tmp_cache = TmpCache(tmp, nothing, nothing, nothing, nothing, nothing)
    return VelocityVerletCache(u, uprev, tmp_cache, k, fsalfirst, half)
end

function alg_cache(
        alg::VelocityVerlet, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return VelocityVerletConstantCache(uEltypeNoUnits(1 // 2))
end

@cache struct LeapfrogDriftKickDriftCache{uType, rateType, uEltypeNoUnits, TmpC <: TmpCache} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    # `TmpCache{*, Nothing, Nothing, Nothing}` — see SymplecticEulerCache.
    tmp_cache::TmpC
    k::rateType
    fsalfirst::rateType
    half::uEltypeNoUnits
end

struct LeapfrogDriftKickDriftConstantCache{uEltypeNoUnits} <: HamiltonConstantCache
    half::uEltypeNoUnits
end

function alg_cache(
        alg::LeapfrogDriftKickDrift, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tmp = zero(rate_prototype)
    k = zero(rate_prototype)
    fsalfirst = zero(rate_prototype)
    half = uEltypeNoUnits(1 // 2)
    # Migrated `tmp` only (was inline); array count unchanged.
    tmp_cache = TmpCache(tmp, nothing, nothing, nothing, nothing, nothing)
    return LeapfrogDriftKickDriftCache(u, uprev, tmp_cache, k, fsalfirst, half)
end

function alg_cache(
        alg::LeapfrogDriftKickDrift, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return LeapfrogDriftKickDriftConstantCache(uEltypeNoUnits(1 // 2))
end

@cache struct VerletLeapfrogCache{uType, rateType, uEltypeNoUnits, TmpC <: TmpCache} <:
    OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    # `TmpCache{*, Nothing, Nothing, Nothing}` — see SymplecticEulerCache.
    tmp_cache::TmpC
    k::rateType
    fsalfirst::rateType
    half::uEltypeNoUnits
end

struct VerletLeapfrogConstantCache{uEltypeNoUnits} <: HamiltonConstantCache
    half::uEltypeNoUnits
end

function alg_cache(
        alg::VerletLeapfrog, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Historically the constructor call swapped the `tmp`/`k` locals, so the
    # `tmp` field held the `zero(rate_prototype)` buffer and `k` held `zero(u)`;
    # that mapping is preserved here. Migrated `tmp` only; array count unchanged.
    tmp = zero(rate_prototype)
    k = zero(u)
    fsalfirst = zero(rate_prototype)
    half = uEltypeNoUnits(1 // 2)
    tmp_cache = TmpCache(tmp, nothing, nothing, nothing, nothing, nothing)
    return VerletLeapfrogCache(u, uprev, tmp_cache, k, fsalfirst, half)
end

function alg_cache(
        alg::VerletLeapfrog, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return VerletLeapfrogConstantCache(uEltypeNoUnits(1 // 2))
end

# ===== Generic symplectic cache for drift-kick methods =====

@cache struct SymplecticGenericCache{uType, rateType, tableauType, TmpC <: TmpCache} <: HamiltonMutableCache
    u::uType
    uprev::uType
    # `TmpCache{*, Nothing, Nothing, Nothing}` — see SymplecticEulerCache.
    tmp_cache::TmpC
    k::rateType
    fsalfirst::rateType
    tab::tableauType
end

# Mapping from algorithm to tableau constructor
_symplectic_tableau(::PseudoVerletLeapfrog, ::Type{T}, ::Type{T2}) where {T, T2} = PseudoVerletLeapfrogConstantCache(T, T2)
_symplectic_tableau(::McAte2, ::Type{T}, ::Type{T2}) where {T, T2} = McAte2ConstantCache(T, T2)
_symplectic_tableau(::Ruth3, ::Type{T}, ::Type{T2}) where {T, T2} = Ruth3ConstantCache(T, T2)
_symplectic_tableau(::McAte3, ::Type{T}, ::Type{T2}) where {T, T2} = McAte3ConstantCache(T, T2)
_symplectic_tableau(::CandyRoz4, ::Type{T}, ::Type{T2}) where {T, T2} = CandyRoz4ConstantCache(T, T2)
_symplectic_tableau(::McAte4, ::Type{T}, ::Type{T2}) where {T, T2} = McAte4ConstantCache(T, T2)
_symplectic_tableau(::CalvoSanz4, ::Type{T}, ::Type{T2}) where {T, T2} = CalvoSanz4ConstantCache(T, T2)
_symplectic_tableau(::McAte42, ::Type{T}, ::Type{T2}) where {T, T2} = McAte42ConstantCache(T, T2)
_symplectic_tableau(::McAte5, ::Type{T}, ::Type{T2}) where {T, T2} = McAte5ConstantCache(T, T2)
_symplectic_tableau(::Yoshida6, ::Type{T}, ::Type{T2}) where {T, T2} = Yoshida6ConstantCache(T, T2)
_symplectic_tableau(::KahanLi6, ::Type{T}, ::Type{T2}) where {T, T2} = KahanLi6ConstantCache(T, T2)
_symplectic_tableau(::McAte8, ::Type{T}, ::Type{T2}) where {T, T2} = McAte8ConstantCache(T, T2)
_symplectic_tableau(::KahanLi8, ::Type{T}, ::Type{T2}) where {T, T2} = KahanLi8ConstantCache(T, T2)
_symplectic_tableau(::SofSpa10, ::Type{T}, ::Type{T2}) where {T, T2} = SofSpa10ConstantCache(T, T2)

const SymplecticGenericAlgorithm = Union{
    PseudoVerletLeapfrog, McAte2, Ruth3, McAte3, CandyRoz4, McAte4,
    CalvoSanz4, McAte42, McAte5, Yoshida6, KahanLi6, McAte8, KahanLi8, SofSpa10,
}

function alg_cache(
        alg::SymplecticGenericAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    # Historically the constructor call swapped the `tmp`/`k` locals, so the
    # `tmp` field held the `zero(rate_prototype)` buffer and `k` held `zero(u)`;
    # that mapping is preserved here. Migrated `tmp` only; array count unchanged.
    tmp = zero(rate_prototype)
    k = zero(u)
    fsalfirst = zero(rate_prototype)
    tab = _symplectic_tableau(
        alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits)
    )
    tmp_cache = TmpCache(tmp, nothing, nothing, nothing, nothing, nothing)
    return SymplecticGenericCache(u, uprev, tmp_cache, k, fsalfirst, tab)
end

function alg_cache(
        alg::SymplecticGenericAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return _symplectic_tableau(
        alg, constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits)
    )
end

function get_fsalfirstlast(
        cache::Union{
            HamiltonMutableCache, VelocityVerletCache, VerletLeapfrogCache,
            SymplecticEulerCache, LeapfrogDriftKickDriftCache,
        }, u
    )
    return (cache.fsalfirst, cache.k)
end
