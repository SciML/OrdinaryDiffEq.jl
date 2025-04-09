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
    P, uType, taylorType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    order::Val{P}
    jet::Function
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

function alg_cache(alg::ExplicitTaylor, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    _, jet_iip = build_jet(f, p, alg.order, length(u))
    utaylor = TaylorDiff.make_seed(u, zero(u), alg.order)
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    ExplicitTaylorCache(alg.order, jet_iip, u, uprev, utaylor, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

get_fsalfirstlast(cache::ExplicitTaylorCache, u) = (cache.u, cache.u)

struct ExplicitTaylorConstantCache{P} <: OrdinaryDiffEqConstantCache
    order::Val{P}
    jet::Function
end
function alg_cache(::ExplicitTaylor{P}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {P, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    if u isa AbstractArray
        jet, _ = build_jet(f, p, Val(P), length(u))
    else
        jet = build_jet(f, p, Val(P))
    end
    ExplicitTaylorConstantCache(Val(P), jet)
end

@cache struct ExplicitTaylorAdaptiveOrderCache{
    uType, taylorType, uNoUnitsType, StageLimiter, StepLimiter,
    Thread} <: OrdinaryDiffEqMutableCache
    min_order::Int
    max_order::Int
    current_order::Ref{Int}
    jets::Vector{Function}
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
    jets = Function[]
    for order in (alg.min_order):(alg.max_order)
        _, jet_iip = build_jet(f, p, Val(order), length(u))
        push!(jets, jet_iip)
    end
    utaylor = TaylorDiff.make_seed(u, zero(u), Val(alg.max_order))
    utilde = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    tmp = zero(u)
    current_order = Ref{Int}(alg.min_order)
    ExplicitTaylorAdaptiveOrderCache(alg.min_order, alg.max_order, current_order,
        jets, u, uprev, utaylor, utilde, tmp, atmp,
        alg.stage_limiter!, alg.step_limiter!, alg.thread)
end

get_fsalfirstlast(cache::ExplicitTaylorAdaptiveOrderCache, u) = (cache.u, cache.u)

struct ExplicitTaylorAdaptiveOrderConstantCache <: OrdinaryDiffEqConstantCache
    min_order::Int
    max_order::Int
    current_order::Ref{Int}
    jets::Vector{Function}
end
function alg_cache(
        alg::ExplicitTaylorAdaptiveOrder, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    jets = Function[]
    for order in (alg.min_order):(alg.max_order)
        if u isa AbstractArray
            jet, _ = build_jet(f, p, Val(order), length(u))
        else
            jet = build_jet(f, p, Val(order))
        end
        push!(jets, jet)
    end
    current_order = Ref{Int}(alg.min_order)
    ExplicitTaylorAdaptiveOrderConstantCache(
        alg.min_order, alg.max_order, current_order, jets)
end
