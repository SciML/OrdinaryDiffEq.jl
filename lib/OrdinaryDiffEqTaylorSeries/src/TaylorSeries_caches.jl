@cache struct ExplicitTaylor2Cache{
        uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
        Thread,
    } <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::rateType
    k2::rateType
    k3::rateType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(
        alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return ExplicitTaylor2Cache(
        u, uprev, k1, k2, k3, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread
    )
end
struct ExplicitTaylor2ConstantCache <: OrdinaryDiffEqConstantCache end
function alg_cache(
        alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return ExplicitTaylor2ConstantCache()
end
# FSAL currently not used, providing dummy implementation to satisfy the interface
get_fsalfirstlast(cache::ExplicitTaylor2Cache, u) = (cache.k1, cache.k1)

@cache struct ExplicitTaylorCache{
        P, jetType, uType, taylorType, coeffType, uNoUnitsType, StageLimiter,
        StepLimiter, Thread,
    } <: OrdinaryDiffEqMutableCache
    order::Val{P}
    jet::jetType
    coeffs::coeffType
    u::uType
    uprev::uType
    utaylor::taylorType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end

function alg_cache(
        alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    _, jet_iip = build_jet(f, p, alg.order, length(u))
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.order)
    coeffs = Vector{eltype(u)}(undef, length(u) * (P + 1))
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    return ExplicitTaylorCache(
        alg.order, jet_iip, coeffs, u, uprev, utaylor, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread
    )
end

get_fsalfirstlast(cache::ExplicitTaylorCache, u) = (nothing, nothing)

struct ExplicitTaylorConstantCache{P, taylorType, uType, tType} <:
    OrdinaryDiffEqConstantCache
    order::Val{P}
    jet::FunctionWrapper{taylorType, Tuple{uType, tType}}
end
function alg_cache(
        alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if u isa AbstractArray
        jet, _ = build_jet(f, p, alg.order, length(u))
    else
        jet = build_jet(f, p, alg.order)
    end
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.order)
    jet_wrapped = FunctionWrapper{typeof(utaylor), Tuple{typeof(u), typeof(t)}}(jet)
    return ExplicitTaylorConstantCache(alg.order, jet_wrapped)
end

mutable struct AdaptiveControllerSnapshot{T}
    q11::T
    errold::T
    dtreject::T
    dtacc::T
    erracc::T
    qold::T
    dt_factor::T
    err1::T
    err2::T
    err3::T
end

function AdaptiveControllerSnapshot(::Type{T}) where {T}
    z = zero(T)
    o = one(T)
    return AdaptiveControllerSnapshot{T}(o, o, o, o, o, o, o, z, z, z)
end

@inline function snapshot_controller!(snap::AdaptiveControllerSnapshot, cache)
    C = typeof(cache)
    hasfield(C, :q11) && (snap.q11 = getfield(cache, :q11))
    hasfield(C, :errold) && (snap.errold = getfield(cache, :errold))
    hasfield(C, :dtreject) && (snap.dtreject = getfield(cache, :dtreject))
    hasfield(C, :dtacc) && (snap.dtacc = getfield(cache, :dtacc))
    hasfield(C, :erracc) && (snap.erracc = getfield(cache, :erracc))
    hasfield(C, :qold) && (snap.qold = getfield(cache, :qold))
    hasfield(C, :dt_factor) && (snap.dt_factor = getfield(cache, :dt_factor))
    if hasfield(C, :err)
        err = getfield(cache, :err)
        snap.err1 = err[1]
        snap.err2 = err[2]
        snap.err3 = err[3]
    end
    return nothing
end

@inline function restore_controller!(cache, snap::AdaptiveControllerSnapshot)
    C = typeof(cache)
    hasfield(C, :q11) && setfield!(cache, :q11, convert(fieldtype(C, :q11), snap.q11))
    hasfield(C, :errold) && setfield!(cache, :errold, convert(fieldtype(C, :errold), snap.errold))
    hasfield(C, :dtreject) &&
        setfield!(cache, :dtreject, convert(fieldtype(C, :dtreject), snap.dtreject))
    hasfield(C, :dtacc) && setfield!(cache, :dtacc, convert(fieldtype(C, :dtacc), snap.dtacc))
    hasfield(C, :erracc) &&
        setfield!(cache, :erracc, convert(fieldtype(C, :erracc), snap.erracc))
    hasfield(C, :qold) && setfield!(cache, :qold, convert(fieldtype(C, :qold), snap.qold))
    hasfield(C, :dt_factor) &&
        setfield!(cache, :dt_factor, convert(fieldtype(C, :dt_factor), snap.dt_factor))
    if hasfield(C, :err)
        err = getfield(cache, :err)
        err[1] = snap.err1
        err[2] = snap.err2
        err[3] = snap.err3
    end
    return nothing
end

@cache struct ExplicitTaylorAdaptiveOrderCache{
        P, Q,
        tType, uType, taylorType, coeffType, uNoUnitsType, StageLimiter, StepLimiter,
        Thread,
    } <: OrdinaryDiffEqMutableCache
    min_order::Val{P}
    max_order::Val{Q}
    current_order::Base.RefValue{Int}
    controller_snapshot::AdaptiveControllerSnapshot{tType}
    jets::Vector{FunctionWrapper{Nothing, Tuple{taylorType, coeffType, uType, tType}}}
    coeffs::Vector{coeffType}
    u::uType
    uprev::uType
    utaylor::taylorType
    utilde::uType
    tmp::uType
    atmp::uNoUnitsType
    stage_limiter!::StageLimiter
    step_limiter!::StepLimiter
    thread::Thread
end
function alg_cache(
        alg::ExplicitTaylorAdaptiveOrder, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_order_value, max_order_value = order_window(alg)
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.max_order)
    coeffType = Vector{eltype(u)}
    jets = FunctionWrapper{
        Nothing, Tuple{typeof(utaylor), coeffType, typeof(u), typeof(t)},
    }[]
    coeffs = coeffType[]
    # every order shares the single `max_order` buffer, so the jets have to fill it
    for order in min_order_value:max_order_value
        jet_iip = build_jet(f, p, Val(order), length(u), alg.max_order)[2]
        push!(jets, jet_iip)
        push!(coeffs, coeffType(undef, length(u) * (order + 1)))
    end
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    current_order = Ref(max_order_value - 1)
    return ExplicitTaylorAdaptiveOrderCache(
        alg.min_order, alg.max_order, current_order, AdaptiveControllerSnapshot(typeof(t)),
        jets, coeffs, u, uprev, utaylor, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread
    )
end

get_fsalfirstlast(cache::ExplicitTaylorAdaptiveOrderCache, u) = (nothing, nothing)

struct ExplicitTaylorAdaptiveOrderConstantCache{
        P, Q, taylorType, uType, tType,
    } <: OrdinaryDiffEqConstantCache
    min_order::Val{P}
    max_order::Val{Q}
    current_order::Base.RefValue{Int}
    controller_snapshot::AdaptiveControllerSnapshot{tType}
    jets::Vector{FunctionWrapper{taylorType, Tuple{uType, tType}}}
end
function alg_cache(
        alg::ExplicitTaylorAdaptiveOrder, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    min_order_value, max_order_value = order_window(alg)
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.max_order) # not used, but needed for type
    jets = FunctionWrapper{typeof(utaylor), Tuple{typeof(u), typeof(t)}}[]
    # the wrapper's return type is pinned at `max_order`, so the jets have to match it
    for order in min_order_value:max_order_value
        if u isa AbstractArray
            jet, _ = build_jet(f, p, Val(order), length(u), alg.max_order)
        else
            jet = build_jet(f, p, Val(order), nothing, alg.max_order)
        end
        push!(jets, jet)
    end
    current_order = Ref(max_order_value - 1)
    return ExplicitTaylorAdaptiveOrderConstantCache(
        alg.min_order, alg.max_order, current_order,
        AdaptiveControllerSnapshot(typeof(t)), jets
    )
end
