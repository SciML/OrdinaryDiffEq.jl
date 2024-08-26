abstract type OrdinaryDiffEqCache <: DiffEqBase.DECache end
abstract type OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache end
abstract type OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache end
struct ODEEmptyCache <: OrdinaryDiffEqConstantCache end
struct ODEChunkCache{CS} <: OrdinaryDiffEqConstantCache end

# Don't worry about the potential alloc on a constant cache
get_fsalfirstlast(cache::OrdinaryDiffEqConstantCache, u) = (zero(u), zero(u))

mutable struct CompositeCache{T, F} <: OrdinaryDiffEqCache
    caches::T
    choice_function::F
    current::Int
end

get_fsalfirstlast(cache::CompositeCache, u) = get_fsalfirstlast(cache.caches[1], u)

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
            args, choice_function, current, u) where {T1, T2, T3, T4, T5, T6, F, uType}
        new{T1, T2, T3, T4, T5, T6, typeof(args), F, uType}(
            args, choice_function, current, u)
    end
end

function get_fsalfirstlast(cache::DefaultCache, u)
    (cache.u, cache.u) # will be overwritten by the cache choice
end

function alg_cache(alg::CompositeAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{V}) where {V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    caches = __alg_cache(alg.algs, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(V))
    CompositeCache(caches, alg.choice_function, 1)
end

function alg_cache(alg::CompositeAlgorithm{CS, Tuple{A1, A2, A3, A4, A5, A6}}, u,
        rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{V}) where {
        CS, V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, A1, A2, A3, A4, A5, A6}
    args = (u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt,
        reltol, p, calck, Val(V))
    argT = map(typeof, args)
    T1 = Base.promote_op(alg_cache, A1, argT...)
    T2 = Base.promote_op(alg_cache, A2, argT...)
    T3 = Base.promote_op(alg_cache, A3, argT...)
    T4 = Base.promote_op(alg_cache, A4, argT...)
    T5 = Base.promote_op(alg_cache, A5, argT...)
    T6 = Base.promote_op(alg_cache, A6, argT...)
    cache = DefaultCache{T1, T2, T3, T4, T5, T6, typeof(alg.choice_function), typeof(u)}(
        args, alg.choice_function, 1, u)
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
    cache
end

# map + closure approach doesn't infer
@generated function __alg_cache(algs::T, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev,
        uprev2, f, t, dt, reltol, p, calck,
        ::Val{V}) where {T <: Tuple, V, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits}
    return Expr(:tuple,
        map(1:length(T.types)) do i
            :(alg_cache(algs[$i], u, rate_prototype, uEltypeNoUnits,
                uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t, dt,
                reltol, p, calck, Val($V)))
        end...)
end

alg_cache(alg::OrdinaryDiffEqAlgorithm, prob, callback::F) where {F} = ODEEmptyCache()

get_chunksize(cache::DiffEqBase.DECache) = error("This cache does not have a chunksize.")
get_chunksize(cache::ODEChunkCache{CS}) where {CS} = CS
