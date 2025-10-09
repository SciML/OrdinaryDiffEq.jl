@cache struct ExplicitTaylor2Cache{
    uType, rateType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
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

function alg_cache(alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    k3 = zero(rate_prototype)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    ExplicitTaylor2Cache(u, uprev, k1, k2, k3, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end
struct ExplicitTaylor2ConstantCache <: OrdinaryDiffEqConstantCache end
function alg_cache(alg::ExplicitTaylor2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    ExplicitTaylor2ConstantCache()
end
# FSAL currently not used, providing dummy implementation to satisfy the interface
get_fsalfirstlast(cache::ExplicitTaylor2Cache, u) = (cache.k1, cache.k1)

@cache struct ExplicitTaylorCache{
    P, tType, uType, taylorType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    order::Val{P}
    jet::FunctionWrapper{Nothing, Tuple{taylorType, uType, tType}}
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

function alg_cache(alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    _, jet_iip = build_jet(f, p, P, length(u))
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.order)
    jet_wrapped = FunctionWrapper{Nothing, Tuple{typeof(utaylor), typeof(u), typeof(t)}}(jet_iip)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    ExplicitTaylorCache(alg.order, jet_wrapped, u, uprev, utaylor, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

get_fsalfirstlast(cache::ExplicitTaylorCache, u) = (cache.u, cache.u)

struct ExplicitTaylorConstantCache{P, taylorType, uType, tType} <:
       OrdinaryDiffEqConstantCache
    order::Val{P}
    jet::FunctionWrapper{taylorType, Tuple{uType, tType}}
end
function alg_cache(alg::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if u isa AbstractArray
        jet, _ = build_jet(f, p, P, length(u))
    else
        jet = build_jet(f, p, P)
    end
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.order) # not used, but needed for type
    jet_wrapped = FunctionWrapper{typeof(utaylor), Tuple{typeof(u), typeof(t)}}(jet)
    ExplicitTaylorConstantCache(alg.order, jet_wrapped)
end

@cache struct ExplicitTaylorAdaptiveOrderCache{P, Q,
    tType, uType, taylorType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    min_order::Val{P}
    max_order::Val{Q}
    current_order::Base.RefValue{Int}
    order_history::Vector{Int}
    jets::Vector{FunctionWrapper{Nothing, Tuple{taylorType, uType, tType}}}
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
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.max_order)
    jets = FunctionWrapper{Nothing, Tuple{typeof(utaylor), typeof(u), typeof(t)}}[]
    min_order_value = get_value(alg.min_order)
    max_order_value = get_value(alg.max_order)
    for order in min_order_value:max_order_value
        jet_iip = build_jet(f, p, order, length(u))[2]
        push!(jets, jet_iip)
    end
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    current_order = Ref(max_order_value - 1)
    order_history = Vector{Int}()
    ExplicitTaylorAdaptiveOrderCache(alg.min_order, alg.max_order, current_order, order_history,
        jets, u, uprev, utaylor, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

get_fsalfirstlast(cache::ExplicitTaylorAdaptiveOrderCache, u) = (cache.u, cache.u)

struct ExplicitTaylorAdaptiveOrderConstantCache{P, Q, taylorType, uType, tType} <:
       OrdinaryDiffEqConstantCache
    min_order::Val{P}
    max_order::Val{Q}
    current_order::Base.RefValue{Int}
    jets::Vector{FunctionWrapper{taylorType, Tuple{uType, tType}}}
end
function alg_cache(
        alg::ExplicitTaylorAdaptiveOrder, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.max_order) # not used, but needed for type
    jets = FunctionWrapper{typeof(utaylor), Tuple{typeof(u), typeof(t)}}[]
    min_order_value = get_value(alg.min_order)
    max_order_value = get_value(alg.max_order)
    for order in min_order_value:max_order_value
        if u isa AbstractArray
            jet, _ = build_jet(f, p, order, length(u))
        else
            jet = build_jet(f, p, order)
        end
        push!(jets, jet)
    end
    current_order = Ref(min_order_value)
    ExplicitTaylorAdaptiveOrderConstantCache(
        alg.min_order, alg.max_order, current_order, jets)
end
