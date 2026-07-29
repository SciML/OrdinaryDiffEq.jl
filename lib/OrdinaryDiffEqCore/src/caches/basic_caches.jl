"""
    OrdinaryDiffEqCache <: SciMLBase.DECache

Abstract supertype of every solver cache. A cache holds the per-solve scratch
state (stage values, temporaries, tableau, nonlinear solver, …) associated with
an algorithm. Concrete caches are built by [`alg_cache`](@ref) and subtype either
[`OrdinaryDiffEqConstantCache`](@ref) or [`OrdinaryDiffEqMutableCache`](@ref).
"""
abstract type OrdinaryDiffEqCache <: SciMLBase.DECache end
"""
    OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache

Abstract supertype for out-of-place ("constant") caches used when the state is
immutable (e.g. `Number`s, `StaticArray`s). These allocate fresh values each step
instead of mutating in place. `is_constant_cache` returns `true` for them.
"""
abstract type OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache end
"""
    OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache

Abstract supertype for in-place ("mutable") caches used when the state is a
mutable array. Their scratch fields are preallocated and reused each step.
"""
abstract type OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache end
struct ODEEmptyCache <: OrdinaryDiffEqConstantCache end
struct ODEChunkCache{CS} <: OrdinaryDiffEqConstantCache end

ismutablecache(cache::OrdinaryDiffEqMutableCache) = true
ismutablecache(cache::OrdinaryDiffEqConstantCache) = false

# Don't worry about the potential alloc on a constant cache
"""
    get_fsalfirstlast(cache, u)

Return the `(fsalfirst, fsallast)` derivative buffers for `cache`, allocating
zeros of the shape of `u` for constant caches. Used to set up FSAL storage when
initializing the integrator.
"""
get_fsalfirstlast(cache::OrdinaryDiffEqConstantCache, u) = (zero(u), zero(u))

"""
    CompositeCache(caches, choice_function, current) <: OrdinaryDiffEqCache

Cache used by [`CompositeAlgorithm`](@ref). Holds the tuple of sub-`caches`, the
`choice_function` that selects the active algorithm, and the index `current` of
the currently-active sub-cache.
"""
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

"""
    DefaultCache <: OrdinaryDiffEqCache

Specialized composite cache used by the automatic default solver. It lazily holds
up to six candidate sub-caches (`cache1`…`cache6`) constructed on demand from
`args`, together with the `choice_function` and the `current` selected index,
avoiding compilation of every candidate up front.
"""
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

"""
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{iip}, verbose)

Construct and return the internal solver cache for `alg`.

This is a developer extension point for solver packages. Each solver sublibrary
defines a method for its algorithm type that allocates stage buffers, tableau
data, nonlinear-solver state, and other per-solve scratch storage. The `Val{iip}`
argument selects the in-place or out-of-place cache branch.

# Examples

Solver packages extend `alg_cache` for their algorithm type:

```julia
import OrdinaryDiffEqCore: alg_cache

function alg_cache(alg::MyAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck, ::Val{iip}, verbose) where {
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, iip}
    return MyAlgorithmCache(u, rate_prototype)
end
```
"""
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
