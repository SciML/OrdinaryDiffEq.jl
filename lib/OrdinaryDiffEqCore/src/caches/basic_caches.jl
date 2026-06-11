abstract type OrdinaryDiffEqCache <: SciMLBase.DECache end
abstract type OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache end
abstract type OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache end
struct ODEEmptyCache <: OrdinaryDiffEqConstantCache end
struct ODEChunkCache{CS} <: OrdinaryDiffEqConstantCache end

ismutablecache(cache::OrdinaryDiffEqMutableCache) = true
ismutablecache(cache::OrdinaryDiffEqConstantCache) = false

# =============================================================================
# Unified scratch buffers for mutable caches.
#
# Historically every mutable cache declared its own ad-hoc scratch fields
# (`tmp`, `utilde`, `atmp`, `linsolve_tmp`, ...) with inconsistent names. That
# made it impossible for shared code (e.g. the initial-dt estimator) to reuse a
# cache's scratch generically. `TmpCache` consolidates the common buffers under
# one struct with stable names, so any cache that carries a `tmp_cache` field
# exposes the same scratch surface.
#
# Parameterized only on the three buffer-array types every cache already has
# (`uType`, `rateType`, `uNoUnitsType`), so adopting it adds no new cache type
# parameter. Fields are concrete (no `Union{T,Nothing}`). Any slot can be opted
# out by parameterizing that slot's type as `Nothing`, which makes the field the
# `nothing` singleton — no array allocated, and `=== nothing` checks fold away
# at compile time:
#   * `TmpCache{uType, rateType, Nothing}` — no unit-less error scratch (`atmp`)
#   * `TmpCache{uType, Nothing, uNoUnitsType}` — no rate scratch, e.g. an
#     explicit method like Tsit5 that never forms a rate-typed linsolve RHS.
#
# NOTE (demo): this is the hand-written core of the larger TmpCache change. Only
# `Tsit5Cache` is wired to it here as a proof of concept; the full migration of
# every cache is intentionally omitted to keep the diff reviewable.
struct TmpCache{uType, rateType, uNoUnitsType}
    tmp::uType          # primary state scratch
    tmp2::uType         # secondary state scratch (e.g. embedded solution / `utilde`)
    atmp::uNoUnitsType  # unit-less state scratch (error norms); `Nothing` if unused
    rate_tmp::rateType  # primary rate scratch (e.g. linsolve RHS / `linsolve_tmp`); `Nothing` if unused
    rate_tmp2::rateType # secondary rate scratch; `Nothing` if unused
end

"""
    build_tmp_cache(u, rate_prototype, uEltypeNoUnits, need_rates::Val = Val(true))

Construct a [`TmpCache`](@ref). The state scratch (`tmp`, `tmp2`) is always
allocated. Two slots are opt-out, controlled at the type level so the result
stays type-stable:

  * pass `Nothing` for `uEltypeNoUnits` to skip the unit-less `atmp` buffer;
  * pass `Val(false)` for `need_rates` to skip both rate buffers (`rate_tmp`,
    `rate_tmp2`) — appropriate for explicit methods that never form a
    rate-typed linsolve RHS.

A skipped slot holds the `nothing` singleton and allocates no array.
"""
function build_tmp_cache(u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Val{need_rates} = Val(true)) where {uEltypeNoUnits, need_rates}
    tmp = zero(u)
    tmp2 = zero(u)
    # `need_rates` is a type parameter, so these ternaries fold at compile time
    # and each specialization returns one concrete `TmpCache` type.
    rate_tmp = need_rates ? zero(rate_prototype) : nothing
    rate_tmp2 = need_rates ? zero(rate_prototype) : nothing
    if uEltypeNoUnits === Nothing
        return TmpCache(tmp, tmp2, nothing, rate_tmp, rate_tmp2)
    else
        atmp = similar(u, uEltypeNoUnits)
        recursivefill!(atmp, false)
        return TmpCache(tmp, tmp2, atmp, rate_tmp, rate_tmp2)
    end
end

# Don't worry about the potential alloc on a constant cache
get_fsalfirstlast(cache::OrdinaryDiffEqConstantCache, u) = (zero(u), zero(u))

mutable struct CompositeCache{T, F} <: OrdinaryDiffEqCache
    caches::T
    choice_function::F
    current::Int
end

function ismutablecache(cache::CompositeCache{T, F}) where {T, F}
    return eltype(T) <: OrdinaryDiffEqMutableCache
end

function get_fsalfirstlast(cache::CompositeCache, u)
    _x = get_fsalfirstlast(cache.caches[1], u)
    if first(_x) !== nothing
        return _x
    else
        return get_fsalfirstlast(cache.caches[2], u)
    end
end

mutable struct DefaultCache{T1, T2, T3, T4, T5, T6, A, F, uType} <: OrdinaryDiffEqCache
    args::A
    choice_function::F
    current::Int
    u::uType
    cache1::T1
    cache2::T2
    cache3::T3
    cache4::T4
    cache5::T5
    cache6::T6
    function DefaultCache{T1, T2, T3, T4, T5, T6, F, uType}(
            args, choice_function, current, u
        ) where {T1, T2, T3, T4, T5, T6, F, uType}
        return new{T1, T2, T3, T4, T5, T6, typeof(args), F, uType}(
            args, choice_function, current, u
        )
    end
end

function get_fsalfirstlast(cache::DefaultCache, u)
    return (cache.u, cache.u) # will be overwritten by the cache choice
end

function ismutablecache(
        cache::DefaultCache{
            T1, T2, T3, T4, T5, T6, A, F, uType,
        }
    ) where {T1, T2, T3, T4, T5, T6, A, F, uType}
    return T1 <: OrdinaryDiffEqMutableCache
end

function alg_cache(
        alg::CompositeAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{V}, verbose
    ) where {V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    caches = __alg_cache(
        alg.algs, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(V), verbose
    )
    return CompositeCache(caches, alg.choice_function, 1)
end

function alg_cache(
        alg::CompositeAlgorithm{CS, Tuple{A1, A2, A3, A4, A5, A6}}, u,
        rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{V}, verbose
    ) where {
        CS, V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, A1, A2, A3, A4, A5, A6,
    }
    args = (
        u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt,
        reltol, p, calck, Val(V), verbose,
    )
    # Core.Typeof to turn uEltypeNoUnits into Type{uEltypeNoUnits} rather than DataType
    argT = map(Core.Typeof, args)
    T1 = Base.promote_op(alg_cache, A1, argT...)
    T2 = Base.promote_op(alg_cache, A2, argT...)
    T3 = Base.promote_op(alg_cache, A3, argT...)
    T4 = Base.promote_op(alg_cache, A4, argT...)
    T5 = Base.promote_op(alg_cache, A5, argT...)
    T6 = Base.promote_op(alg_cache, A6, argT...)
    cache = DefaultCache{T1, T2, T3, T4, T5, T6, typeof(alg.choice_function), typeof(u)}(
        args, alg.choice_function, 1, u
    )
    algs = alg.algs
    # If the type is a bitstype we need to initialize it correctly here since isdefined will always return true.
    if isbitstype(T1)
        cache.cache1 = alg_cache(algs[1], args...)
    end
    if isbitstype(T2)
        cache.cache2 = alg_cache(algs[2], args...)
    end
    if isbitstype(T3)
        cache.cache3 = alg_cache(algs[3], args...)
    end
    if isbitstype(T4)
        cache.cache4 = alg_cache(algs[4], args...)
    end
    if isbitstype(T5)
        cache.cache5 = alg_cache(algs[5], args...)
    end
    if isbitstype(T6)
        cache.cache6 = alg_cache(algs[6], args...)
    end
    return cache
end

# map + closure approach doesn't infer
@generated function __alg_cache(
        algs::T, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
        uprev2, f, t, dt, reltol, p, calck,
        ::Val{V}, verbose
    ) where {
        T <: Tuple, V, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits,
    }
    return Expr(
        :tuple,
        map(1:length(T.types)) do i
            :(
                alg_cache(
                    algs[$i], u, rate_prototype, uEltypeNoUnits,
                    uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt,
                    reltol, p, calck, Val($V), verbose
                )
            )
        end...
    )
end

alg_cache(alg::OrdinaryDiffEqAlgorithm, prob, callback::F) where {F} = ODEEmptyCache()

get_chunksize(cache::SciMLBase.DECache) = error("This cache does not have a chunksize.")
get_chunksize(cache::ODEChunkCache{CS}) where {CS} = CS
