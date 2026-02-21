abstract type DAEBDFMutableCache <: OrdinaryDiffEqMutableCache end
function get_fsalfirstlast(cache::DAEBDFMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

@cache mutable struct DImplicitEulerCache{uType, rateType, uNoUnitsType, N} <:
    DAEBDFMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    atmp::uNoUnitsType
    k₁::rateType
    k₂::rateType
    nlsolver::N
end

# Not FSAL
get_fsalfirstlast(cache::DImplicitEulerCache, u) = (nothing, nothing)

mutable struct DImplicitEulerConstantCache{N} <: OrdinaryDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::DImplicitEuler, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    α = 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(false), verbose
    )

    return DImplicitEulerConstantCache(nlsolver)
end

function alg_cache(
        alg::DImplicitEuler, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    α = 1
    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(true), verbose
    )

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    return DImplicitEulerCache(u, uprev, uprev2, atmp, k₁, k₂, nlsolver)
end

@cache mutable struct DABDF2ConstantCache{N, dtType, rate_prototype} <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    eulercache::DImplicitEulerConstantCache
    dtₙ₋₁::dtType
    fsalfirstprev::rate_prototype
end

function alg_cache(
        alg::DABDF2, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(1) // 1, 1
    α = Int64(1) // 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(false), verbose
    )
    eulercache = DImplicitEulerConstantCache(nlsolver)

    dtₙ₋₁ = one(dt)
    fsalfirstprev = rate_prototype

    return DABDF2ConstantCache(nlsolver, eulercache, dtₙ₋₁, fsalfirstprev)
end

@cache mutable struct DABDF2Cache{uType, rateType, uNoUnitsType, N, dtType} <:
    DAEBDFMutableCache
    uₙ::uType
    uₙ₋₁::uType
    uₙ₋₂::uType
    fsalfirst::rateType
    fsalfirstprev::rateType
    atmp::uNoUnitsType
    nlsolver::N
    eulercache::DImplicitEulerCache
    dtₙ₋₁::dtType
end

function alg_cache(
        alg::DABDF2, du, u, res_prototype, rate_prototype,
        ::Type{uEltypeNoUnits}, ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(1) // 1, 1
    α = Int64(1) // 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, res_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, α, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    fsalfirstprev = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    k₁ = zero(rate_prototype)
    k₂ = zero(rate_prototype)

    eulercache = DImplicitEulerCache(u, uprev, uprev2, atmp, k₁, k₂, nlsolver)

    dtₙ₋₁ = one(dt)

    return DABDF2Cache(
        u, uprev, uprev2, fsalfirst, fsalfirstprev, atmp,
        nlsolver, eulercache, dtₙ₋₁
    )
end

@cache mutable struct DFBDFConstantCache{
        MO, N, tsType, tType, uType, uuType, coeffType,
        EEstType, rType, wType,
    } <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    ts::tsType
    ts_tmp::tsType
    t_old::tType
    u_history::uuType
    order::Int
    prev_order::Int
    u₀::uType
    u_corrector::uuType
    bdf_coeffs::coeffType
    max_order::Val{MO}
    nconsteps::Int
    consfailcnt::Int
    terkm2::EEstType
    terkm1::EEstType
    terk::EEstType
    terkp1::EEstType
    r::rType
    weights::wType
    iters_from_event::Int
end

function alg_cache(
        alg::DFBDF{MO}, du, u, res_prototype, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits,
        uprev, uprev2, f, t, dt, reltol, p, calck, ::Val{false}, verbose
    ) where {MO}
    γ, c = 1.0, 1.0
    max_order = MO
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    bdf_coeffs = SA[
        1 -1 0 0 0 0;
        Int64(2) // 3 -Int64(4) // 3 Int64(1) // 3 0 0 0;
        Int64(6) // 11 -Int64(18) // 11 Int64(9) // 11 -Int64(2) // 11 0 0;
        Int64(12) // 25 -Int64(48) // 25 Int64(36) // 25 -Int64(16) // 25 Int64(3) // 25 0;
        Int64(60) // 137 -Int64(300) // 137 Int64(300) // 137 -Int64(200) // 137 Int64(75) // 137 -Int64(12) // 137
    ]
    ts = zero(Vector{typeof(t)}(undef, max_order + 2)) #ts is the successful past points, it will be updated after successful step
    ts_tmp = similar(ts)

    u_history = zero(Matrix{eltype(u)}(undef, length(u), max_order + 2))
    order = 1
    prev_order = 1
    u_corrector = similar(u_history)
    recursivefill!(u_corrector, zero(eltype(u)))
    recursivefill!(u_history, zero(eltype(u_history)))
    terkm2 = tTypeNoUnits(1)
    terkm1 = tTypeNoUnits(1)
    terk = tTypeNoUnits(1)
    terkp1 = tTypeNoUnits(1)
    r = zero(Vector{typeof(t)}(undef, max_order + 2))
    weights = zero(Vector{typeof(t)}(undef, max_order + 2))
    weights[1] = 1
    nconsteps = 0
    consfailcnt = 0
    t_old = zero(t)
    iters_from_event = 0
    u₀ = zero(u)

    return DFBDFConstantCache(
        nlsolver, ts, ts_tmp, t_old, u_history, order, prev_order, u₀,
        u_corrector, bdf_coeffs, Val(MO), nconsteps, consfailcnt, terkm2,
        terkm1, terk, terkp1, r, weights, iters_from_event
    )
end

@cache mutable struct DFBDFCache{
        MO, N, rateType, uNoUnitsType, tsType, tType, uType,
        uuType, coeffType, EEstType, rType, wType,
    } <:
    DAEBDFMutableCache
    fsalfirst::rateType
    nlsolver::N
    ts::tsType
    ts_tmp::tsType
    t_old::tType
    u_history::uuType
    order::Int
    prev_order::Int
    u_corrector::uuType
    u₀::uType
    bdf_coeffs::coeffType
    max_order::Val{MO}
    nconsteps::Int
    consfailcnt::Int
    tmp::uType
    atmp::uNoUnitsType
    terkm2::EEstType
    terkm1::EEstType
    terk::EEstType
    terkp1::EEstType
    terk_tmp::uType
    terkp1_tmp::uType
    r::rType
    weights::wType
    equi_ts::tsType
    iters_from_event::Int
    dense::Vector{uType}
end

function alg_cache(
        alg::DFBDF{MO}, du, u, res_prototype, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {MO}
    γ, c = 1.0, 1.0
    fsalfirst = zero(rate_prototype)
    max_order = MO
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    #=bdf_coeffs = SA[1 -1 0 0 0 0 ;
                    3//2 -2 1//2 0 0 0 ;
                    11//6 -3 3//2 -1//3  0 0 ;
                    25//12 -4 3 -4//3 1//4 0 ;
                    137//60 -5 5 -10//3 5//4 -1//5]=#
    bdf_coeffs = SA[
        1 -1 0 0 0 0;
        Int64(2) // 3 -Int64(4) // 3 Int64(1) // 3 0 0 0;
        Int64(6) // 11 -Int64(18) // 11 Int64(9) // 11 -Int64(2) // 11 0 0;
        Int64(12) // 25 -Int64(48) // 25 Int64(36) // 25 -Int64(16) // 25 Int64(3) // 25 0;
        Int64(60) // 137 -Int64(300) // 137 Int64(300) // 137 -Int64(200) // 137 Int64(75) // 137 -Int64(12) // 137
    ]
    ts = Vector{typeof(t)}(undef, max_order + 2) #ts is the successful past points, it will be updated after successful step
    u_history = Matrix{eltype(u)}(undef, length(u), max_order + 2)
    order = 1
    prev_order = 1
    u_corrector = similar(u_history)
    recursivefill!(ts, zero(t))
    recursivefill!(u_corrector, zero(eltype(u)))
    recursivefill!(u_history, zero(eltype(u_history)))
    terkm2 = tTypeNoUnits(1)
    terkm1 = tTypeNoUnits(1)
    terk = tTypeNoUnits(1)
    terkp1 = tTypeNoUnits(1)
    terk_tmp = similar(u)
    terkp1_tmp = similar(u)
    r = Vector{typeof(t)}(undef, max_order + 2)
    weights = Vector{typeof(t)}(undef, max_order + 2)
    recursivefill!(r, zero(t))
    recursivefill!(weights, zero(t))
    weights[1] = 1
    nconsteps = 0
    consfailcnt = 0
    t_old = zero(t)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, zero(uEltypeNoUnits))
    u₀ = similar(u)
    equi_ts = similar(ts)
    tmp = similar(u)
    ts_tmp = similar(ts)
    iters_from_event = 0

    dense = [zero(u) for _ in 1:(2 * (max_order + 1))]

    return DFBDFCache(
        fsalfirst, nlsolver, ts, ts_tmp, t_old, u_history, order, prev_order,
        u_corrector, u₀, bdf_coeffs, Val(MO), nconsteps, consfailcnt, tmp, atmp,
        terkm2, terkm1, terk, terkp1, terk_tmp, terkp1_tmp, r, weights, equi_ts,
        iters_from_event, dense
    )
end
