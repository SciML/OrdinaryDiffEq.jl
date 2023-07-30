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

if isdefined(Base, :Experimental) && isdefined(Base.Experimental, :silence!)
    Base.Experimental.silence!(CompositeCache)
end

function alg_cache(alg::CompositeAlgorithm{Tuple{T1, T2}, F}, u, rate_prototype,
    ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
    ::Type{tTypeNoUnits}, uprev,
    uprev2, f, t, dt, reltol, p, calck,
    ::Val{V}) where {T1, T2, F, V, uEltypeNoUnits, uBottomEltypeNoUnits,
    tTypeNoUnits}
    caches = (alg_cache(alg.algs[1], u, rate_prototype, uEltypeNoUnits,
            uBottomEltypeNoUnits,
            tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(V)),
        alg_cache(alg.algs[2], u, rate_prototype, uEltypeNoUnits,
            uBottomEltypeNoUnits,
            tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(V)))
    CompositeCache(caches, alg.choice_function, 1)
end

function alg_cache(alg::CompositeAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
    ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
    dt, reltol, p, calck,
    ::Val{V}) where {V, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    caches = __alg_cache(alg.algs, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(V))
    CompositeCache(caches, alg.choice_function, 1)
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
