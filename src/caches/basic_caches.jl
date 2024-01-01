abstract type OrdinaryDiffEqCache <: DiffEqBase.DECache end
abstract type OrdinaryDiffEqConstantCache <: OrdinaryDiffEqCache end
abstract type OrdinaryDiffEqMutableCache <: OrdinaryDiffEqCache end
struct ODEEmptyCache <: OrdinaryDiffEqConstantCache end
struct ODEChunkCache{CS} <: OrdinaryDiffEqConstantCache end

mutable struct CompositeCache{T, F} <: OrdinaryDiffEqCache
    caches::T
    choice_function::F
    current::Int
end

TruncatedStacktraces.@truncate_stacktrace CompositeCache 1

mutable struct DefaultCache{T1, T2, T3, T4, T5, T6, A, F} <: OrdinaryDiffEqCache
    args::A
    choice_function::F
    current::Int
    cache1::T1
    cache2::T2
    cache3::T3
    cache4::T4
    cache5::T5
    cache6::T6
    function DefaultCache{T1, T2, T3, T4, T5, T6, F}(args, choice_function, current) where {T1, T2, T3, T4, T5, T6, F}
        new{T1, T2, T3, T4, T5, T6, typeof(args), F}(args, choice_function, current)
    end
end

TruncatedStacktraces.@truncate_stacktrace DefaultCache 1

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, :silence!)
    Base.Experimental.silence!(CompositeCache)
    Base.Experimental.silence!(DefaultCache)
end

function alg_cache(alg::CompositeAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{V}) where {V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    caches = __alg_cache(alg.algs, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(V))
    CompositeCache{typeof(caches), typeof(alg.choice_function)}(
            caches, alg.choice_function, 1)
end

function alg_cache(alg::CompositeAlgorithm{Tuple{A1, A2, A3, A4, A5, A6}}, u,
    rate_prototype, ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
    uprev, uprev2, f, t, dt, reltol, p, calck,
    ::Val{V}) where {V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, A1, A2, A3, A4, A5, A6}

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
    DefaultCache{T1, T2, T3, T4, T5, T6, typeof(alg.choice_function)}(args, alg.choice_function, 1)
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

@cache struct FunctionMapCache{uType, rateType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::rateType
end

function alg_cache(alg::FunctionMap, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    FunctionMapCache(u, uprev,
        FunctionMap_scale_by_time(alg) ? rate_prototype :
        (eltype(u) <: Enum ? copy(u) : zero(u)))
end

struct FunctionMapConstantCache <: OrdinaryDiffEqConstantCache end

function alg_cache(alg::FunctionMap, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    FunctionMapConstantCache()
end

@cache struct ExplicitRKCache{uType, rateType, uNoUnitsType, TabType} <:
              OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    utilde::rateType
    atmp::uNoUnitsType
    fsalfirst::rateType
    fsallast::rateType
    kk::Vector{rateType}
    tab::TabType
end

TruncatedStacktraces.@truncate_stacktrace ExplicitRKCache 1

function alg_cache(alg::ExplicitRK, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    kk = Vector{typeof(rate_prototype)}(undef, 0)
    for i in 1:(alg.tableau.stages)
        push!(kk, zero(rate_prototype))
    end
    fsalfirst = kk[1]
    if isfsal(alg.tableau)
        fsallast = kk[end]
    else
        fsallast = zero(rate_prototype)
    end
    utilde = zero(rate_prototype)
    tmp = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tab = ExplicitRKConstantCache(alg.tableau, rate_prototype)
    ExplicitRKCache(u, uprev, tmp, utilde, atmp, fsalfirst, fsallast, kk, tab)
end

struct ExplicitRKConstantCache{MType, VType, KType} <: OrdinaryDiffEqConstantCache
    A::MType
    c::VType
    α::VType
    αEEst::VType
    stages::Int
    kk::KType
end

function ExplicitRKConstantCache(tableau, rate_prototype)
    @unpack A, c, α, αEEst, stages = tableau
    A = copy(A') # Transpose A to column major looping
    kk = Array{typeof(rate_prototype)}(undef, stages) # Not ks since that's for integrator.opts.dense
    αEEst = isempty(αEEst) ? αEEst : α .- αEEst
    ExplicitRKConstantCache(A, c, α, αEEst, stages, kk)
end

function alg_cache(alg::ExplicitRK, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ExplicitRKConstantCache(alg.tableau, rate_prototype)
end

get_chunksize(cache::DiffEqBase.DECache) = error("This cache does not have a chunksize.")
get_chunksize(cache::ODEChunkCache{CS}) where {CS} = CS
