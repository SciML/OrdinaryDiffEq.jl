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
# parameter. Fields are concrete (no `Union{T,Nothing}`).
#
# The slots are deliberately sized to cover the initial-dt estimator's
# scratch needs (see `_ode_initdt_iip`): unit-less scale buffers, two state
# temporaries, and two rate temporaries. Slots are populated ALIASING-FIRST to
# keep the total cache size capped (many methods are sensitive to the exact
# number of u-sized vectors they allocate — low-storage RK in particular):
#
#   1. Migration: the cache's existing scratch fields (`tmp`, `utilde`, `atmp`,
#      `weight`) move into the slots — the field is deleted, single name, zero
#      new arrays.
#   2. Donor aliasing: a slot with no migrated field may alias another buffer
#      the cache already owns, PROVIDED the donor is dead between steps AND not
#      read by dense output / `_ode_addsteps!` (initdt may scribble on slots
#      mid-solve via callback-triggered `auto_dt_reset!`; interpolation stage
#      arrays are therefore never legal donors). Construct via the raw
#      `TmpCache(...)` constructor in that case (see OrdinaryDiffEqNordsieck
#      for the pattern).
#   3. Fresh allocation: only when the user opts in via the algorithm's
#      `preallocate_initdt_buffers` option (default `false`); otherwise the
#      slot is `nothing` and `initdt` falls back to allocating at call time.
#
# A slot is opted out by parameterizing its type as `Nothing`: the field becomes
# the `nothing` singleton, no array is allocated, and `=== nothing` checks fold
# away at compile time.
#   * `TmpCache{uType, Nothing, rateType, Nothing, Nothing}` — no secondary
#     state (`tmp2`) and no unit-less scratch (`atmp`/`weight`); e.g. nlsolver
#     caches that only need a `tmp` and an `atmp`.
#   * `TmpCache{uType, uType, Nothing, uNoUnitsType, Nothing}` — no rate scratch.
#     The rate buffers exist only so `initdt` can run allocation-free; a cache
#     may skip them and let `initdt` fall back to allocating (see the
#     `preallocate_initdt_buffers` algorithm option). `initdt` reuses whatever
#     buffers are present.
#   * `weight` is a second unit-less slot (e.g. Rosenbrock's linear-solve
#     weighting); most caches opt out.
#
# `tmp` and `tmp2` carry SEPARATE type parameters (`uType`, `uType2`) even
# though both are state-typed: a cache that has a `tmp` but no `utilde`-style
# secondary buffer sets `tmp2 = nothing`, which a shared parameter would forbid
# (the diagonal rule rejects `TmpCache(array, nothing, ...)`).
struct TmpCache{uType, uType2, rateType, uNoUnitsType, weightType}
    tmp::uType          # primary state scratch
    tmp2::uType2        # secondary state scratch (e.g. embedded solution / `utilde`); `Nothing` if unused
    atmp::uNoUnitsType  # unit-less state scratch (error norms); `Nothing` if unused
    weight::weightType  # secondary unit-less scratch (e.g. linsolve weights); `Nothing` if unused
    rate_tmp::rateType  # primary rate scratch (initdt `f₁`); `Nothing` if unused
    rate_tmp2::rateType # secondary rate scratch (initdt mass-matrix `ftmp`); `Nothing` if unused
end

# Single source of truth for the resizable buffers a `TmpCache` carries. The
# `@cache` macro splats this into `full_cache` (so `resize!` & friends see the
# sub-buffers), which means adding a slot to the struct only requires updating
# this tuple — the macro never has to change. Keep in sync with the struct
# fields above; opted-out slots appear as `nothing` and are skipped by
# `full_cache` consumers.
function tmp_cache_buffers(tc::TmpCache)
    return (tc.tmp, tc.tmp2, tc.atmp, tc.weight, tc.rate_tmp, tc.rate_tmp2)
end

"""
    build_tmp_cache(u, rate_prototype, uEltypeNoUnits,
        need_rates::Val = Val(false), need_weight::Val = Val(false))

Construct a [`TmpCache`](@ref) for the common cache layout, allocating the
state scratch (`tmp`, `tmp2`) and the unit-less `atmp`. This is the helper for
caches that previously allocated `tmp`/`utilde`/`atmp` inline — using it is
net-zero in array count. A cache whose former scratch layout differs (no
`utilde`, donor-aliased rate buffers, ...) should call the raw
`TmpCache(tmp, tmp2, atmp, weight, rate_tmp, rate_tmp2)` constructor instead,
passing its own arrays (or `nothing`) per slot — see the aliasing policy above
the struct definition.

The optional slots are controlled at the type level so the result is
type-stable:

  * pass `Nothing` for `uEltypeNoUnits` to skip the unit-less buffers (`atmp`
    and `weight`);
  * pass `Val(true)` for `need_rates` to allocate the rate buffers (`rate_tmp`,
    `rate_tmp2`). These exist only to make `initdt` allocation-free; they
    default OFF so the cache's footprint is unchanged unless the user opts in
    via the `preallocate_initdt_buffers` algorithm option — wire it through as
    `Val(preallocate_initdt_buffers(alg))`;
  * pass `Val(true)` for `need_weight` to allocate the secondary unit-less
    `weight` buffer (e.g. Rosenbrock's linear-solve weighting).

A skipped slot holds the `nothing` singleton and allocates no array.
"""
function build_tmp_cache(u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Val{need_rates} = Val(false), ::Val{need_weight} = Val(false)
    ) where {uEltypeNoUnits, need_rates, need_weight}
    tmp = zero(u)
    tmp2 = zero(u)
    # `need_rates`/`need_weight` are type parameters, so these ternaries fold at
    # compile time and each specialization returns a single concrete `TmpCache`.
    rate_tmp = need_rates ? zero(rate_prototype) : nothing
    rate_tmp2 = need_rates ? zero(rate_prototype) : nothing
    if uEltypeNoUnits === Nothing
        return TmpCache(tmp, tmp2, nothing, nothing, rate_tmp, rate_tmp2)
    else
        atmp = similar(u, uEltypeNoUnits)
        recursivefill!(atmp, false)
        weight = if need_weight
            w = similar(u, uEltypeNoUnits)
            recursivefill!(w, false)
            w
        else
            nothing
        end
        return TmpCache(tmp, tmp2, atmp, weight, rate_tmp, rate_tmp2)
    end
end

"""
    preallocate_initdt_buffers(alg)::Bool

Whether `alg`'s cache should allocate dedicated rate-scratch buffers so the
initial-`dt` estimator (and `reinit!`/`auto_dt_reset!`) runs allocation-free
even where no safe donor buffer exists to alias. Reads the algorithm's
`preallocate_initdt_buffers` field when it has one (the user-facing knob), and
defaults to `false` otherwise: the default cache footprint is always identical
to the historical one — extra arrays are strictly opt-in. `hasfield` on a
concrete algorithm type is a compile-time constant, so this folds away.
"""
function preallocate_initdt_buffers(alg)
    if hasfield(typeof(alg), :preallocate_initdt_buffers)
        return alg.preallocate_initdt_buffers
    else
        return false
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

# =============================================================================
# Accessor for the unified scratch struct. Shared code (the initial-dt
# estimator, `reinit!`/`auto_dt_reset!`) reaches a cache's `TmpCache` through
# this rather than through `get_tmp_cache` — the latter keeps its historical
# positional-tuple contract for callbacks/DelayDiffEq. Returns `nothing` for
# constant caches and any cache without a `tmp_cache` field; callers fall back
# to allocating. `hasfield` on a concrete cache type is a compile-time
# constant, so the branch folds away.
@inline initdt_tmp_cache(cache) =
    hasfield(typeof(cache), :tmp_cache) ? cache.tmp_cache : nothing

# Composite/default solves forward to the currently selected method's cache so
# they get allocation-free initdt too. Runtime indexing here is type-unstable,
# but this runs once per `init`/`auto_dt_reset!`, not in the step loop.
initdt_tmp_cache(cache::CompositeCache) = initdt_tmp_cache(cache.caches[cache.current])

function initdt_tmp_cache(cache::DefaultCache)
    # `cache1..cache6` are lazily constructed; guard with `isdefined` since
    # initdt can run before the choice function has instantiated the current
    # cache (then we just fall back to allocating).
    n = cache.current
    if n == 1 && isdefined(cache, :cache1)
        return initdt_tmp_cache(cache.cache1)
    elseif n == 2 && isdefined(cache, :cache2)
        return initdt_tmp_cache(cache.cache2)
    elseif n == 3 && isdefined(cache, :cache3)
        return initdt_tmp_cache(cache.cache3)
    elseif n == 4 && isdefined(cache, :cache4)
        return initdt_tmp_cache(cache.cache4)
    elseif n == 5 && isdefined(cache, :cache5)
        return initdt_tmp_cache(cache.cache5)
    elseif n == 6 && isdefined(cache, :cache6)
        return initdt_tmp_cache(cache.cache6)
    else
        return nothing
    end
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
