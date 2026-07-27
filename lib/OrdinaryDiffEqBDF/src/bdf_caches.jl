abstract type BDFMutableCache <: OrdinaryDiffEqMutableCache end
function get_fsalfirstlast(cache::BDFMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

@cache mutable struct ABDF2ConstantCache{N, dtType, rate_prototype, EC} <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    eulercache::EC
    dtₙ₋₁::dtType
    fsalfirstprev::rate_prototype
end

function alg_cache(
        alg::ABDF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(2) // 3, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    ie_tab = ImplicitEulerESDIRKIMEXTableau(
        constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits)
    )
    eulercache = ESDIRKIMEXConstantCache(nlsolver, ie_tab)

    dtₙ₋₁ = one(dt)
    fsalfirstprev = rate_prototype

    return ABDF2ConstantCache(nlsolver, eulercache, dtₙ₋₁, fsalfirstprev)
end

@cache mutable struct ABDF2Cache{
        uType, rateType, N, dtType, StepLimiter, EC, TmpC <: TmpCache,
    } <: BDFMutableCache
    uₙ::uType
    uₙ₋₁::uType
    uₙ₋₂::uType
    fsalfirst::rateType
    fsalfirstprev::rateType
    zₙ₋₁::uType
    # Unified scratch: only the cache-level `atmp` migrated here (the Newton
    # buffers live on the nlsolver and are off limits). `tmp`/`tmp2`/`weight`
    # are `nothing`; the rate slots stay opted out since ABDF2's positional
    # constructor exposes no `preallocate_initdt_buffers` knob:
    # `TmpCache{Nothing, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    nlsolver::N
    eulercache::EC
    dtₙ₋₁::dtType
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::ABDF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(2) // 3, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    fsalfirstprev = zero(rate_prototype)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    algebraic_vars = f.mass_matrix === I ? nothing :
        find_algebraic_vars_eqs(f.mass_matrix)[1]

    ie_tab = ImplicitEulerESDIRKIMEXTableau(
        constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits)
    )
    zs = [nlsolver.z]
    ks = Vector{Nothing}()
    eulercache = ESDIRKIMEXCache(
        u, uprev, fsalfirst, zs, ks, atmp, nlsolver, ie_tab, alg.step_limiter!,
        uprev2, algebraic_vars, nothing, nothing
    )

    dtₙ₋₁ = one(dt)
    zₙ₋₁ = zero(u)

    # `atmp` (shared with the bootstrap eulercache, preserving the historical
    # aliasing) is the only migrated scratch field; everything else stays
    # `nothing`, so the array count matches the historical cache exactly.
    tmp_cache = TmpCache(nothing, nothing, atmp, nothing, nothing, nothing)

    return ABDF2Cache(
        u, uprev, uprev2, fsalfirst, fsalfirstprev, zₙ₋₁, tmp_cache,
        nlsolver, eulercache, dtₙ₋₁, alg.step_limiter!
    )
end

# SBDF

@cache mutable struct SBDFConstantCache{rateType, N, uType} <: OrdinaryDiffEqConstantCache
    cnt::Int
    ark::Bool
    k2::rateType
    nlsolver::N
    uprev2::uType
    uprev4::uType
    uprev3::uType
    k₁::rateType
    k₂::rateType
    k₃::rateType
    du₁::rateType
    du₂::rateType
end

@cache mutable struct SBDFCache{uType, rateType, N} <: BDFMutableCache
    cnt::Int
    ark::Bool
    u::uType
    uprev::uType
    fsalfirst::rateType
    nlsolver::N
    uprev2::uType
    uprev3::uType
    uprev4::uType
    k₁::rateType
    k₂::rateType
    k₃::rateType
    du₁::rateType
    du₂::rateType
end

function alg_cache(
        alg::SBDF, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(1) // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    k2 = rate_prototype
    k₁ = rate_prototype
    k₂ = rate_prototype
    k₃ = rate_prototype
    du₁ = rate_prototype
    du₂ = rate_prototype

    uprev2 = u
    uprev3 = u
    uprev4 = u

    return SBDFConstantCache(
        1, alg.ark, k2, nlsolver, uprev2, uprev3, uprev4, k₁, k₂, k₃, du₁,
        du₂
    )
end

function alg_cache(
        alg::SBDF, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = Int64(1) // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    order = alg.order

    k₁ = zero(rate_prototype)
    k₂ = order >= 3 ? zero(rate_prototype) : k₁
    k₃ = order == 4 ? zero(rate_prototype) : k₁
    du₁ = zero(rate_prototype)
    du₂ = zero(rate_prototype)

    uprev2 = zero(u)
    uprev3 = order >= 3 ? zero(u) : uprev2
    uprev4 = order == 4 ? zero(u) : uprev2

    return SBDFCache(
        1, alg.ark, u, uprev, fsalfirst, nlsolver, uprev2, uprev3, uprev4, k₁, k₂, k₃,
        du₁, du₂
    )
end

# QNDF1

@cache mutable struct QNDF1ConstantCache{
        N,
        coefType,
        coefType1,
        coefType2,
        dtType,
        uType,
    } <: OrdinaryDiffEqConstantCache
    nlsolver::N
    D::coefType1
    D2::coefType2
    R::coefType
    U::coefType
    uprev2::uType
    dtₙ₋₁::dtType
end

@cache mutable struct QNDF1Cache{
        uType, rateType, coefType, coefType1, coefType2,
        N, dtType, StepLimiter, TmpC <: TmpCache,
    } <: BDFMutableCache
    uprev2::uType
    fsalfirst::rateType
    D::coefType1
    D2::coefType2
    R::coefType
    U::coefType
    # Unified scratch: `utilde` migrated into `tmp2` and `atmp` into `atmp`.
    # The `tmp` slot donor-aliases `nlsolver.tmp` (recomputed write-first every
    # step, dead between steps, never read by dense output — master's initdt
    # already borrows it via `get_tmp_cache` for Newton algorithms). Rate slots
    # stay opted out (no `preallocate_initdt_buffers` knob on QNDF1):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    nlsolver::N
    dtₙ₋₁::dtType
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::QNDF1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = zero(inv((1 - alg.kappa))), 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    uprev2 = u
    dtₙ₋₁ = zero(t)

    D = fill(zero(u), 1, 1)
    D2 = fill(zero(u), 1, 2)
    R = fill(zero(t), 1, 1)
    U = fill(zero(t), 1, 1)

    U!(1, U)

    return QNDF1ConstantCache(nlsolver, D, D2, R, U, uprev2, dtₙ₋₁)
end

function alg_cache(
        alg::QNDF1, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = zero(inv((1 - alg.kappa))), 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    D = Array{typeof(u)}(undef, 1, 1)
    D2 = Array{typeof(u)}(undef, 1, 2)

    R = fill(zero(t), 1, 1)
    U = fill(zero(t), 1, 1)

    D[1] = zero(u)
    D2[1] = zero(u)
    D2[2] = zero(u)

    U!(1, U)

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    utilde = zero(u)
    uprev2 = zero(u)
    dtₙ₋₁ = zero(dt)

    # Migrated `utilde`/`atmp` plus the `nlsolver.tmp` donor alias — zero new
    # arrays vs. the historical cache.
    tmp_cache = TmpCache(nlsolver.tmp, utilde, atmp, nothing, nothing, nothing)

    return QNDF1Cache(
        uprev2, fsalfirst, D, D2, R, U, tmp_cache, nlsolver, dtₙ₋₁, alg.step_limiter!
    )
end

# QNDF2

@cache mutable struct QNDF2ConstantCache{
        N,
        coefType,
        coefType1,
        coefType2,
        uType,
        dtType,
    } <: OrdinaryDiffEqConstantCache
    nlsolver::N
    D::coefType1
    D2::coefType2
    R::coefType
    U::coefType
    uprev2::uType
    uprev3::uType
    dtₙ₋₁::dtType
    dtₙ₋₂::dtType
end

@cache mutable struct QNDF2Cache{
        uType, rateType, coefType, coefType1, coefType2,
        N, dtType, StepLimiter, TmpC <: TmpCache,
    } <: BDFMutableCache
    uprev2::uType
    uprev3::uType
    fsalfirst::rateType
    D::coefType1
    D2::coefType2
    R::coefType
    U::coefType
    # Unified scratch: `utilde` migrated into `tmp2` and `atmp` into `atmp`.
    # The `tmp` slot donor-aliases `nlsolver.tmp` (recomputed write-first every
    # step, dead between steps, never read by dense output — master's initdt
    # already borrows it via `get_tmp_cache` for Newton algorithms). Rate slots
    # stay opted out (no `preallocate_initdt_buffers` knob on QNDF2):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    nlsolver::N
    dtₙ₋₁::dtType
    dtₙ₋₂::dtType
    step_limiter!::StepLimiter
end

function alg_cache(
        alg::QNDF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = zero(inv((1 - alg.kappa))), 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    uprev2 = u
    uprev3 = u
    dtₙ₋₁ = zero(t)
    dtₙ₋₂ = zero(t)

    D = fill(zero(u), 1, 2)
    D2 = fill(zero(u), 1, 3)
    R = fill(zero(t), 2, 2)
    U = fill(zero(t), 2, 2)

    U!(2, U)

    return QNDF2ConstantCache(nlsolver, D, D2, R, U, uprev2, uprev3, dtₙ₋₁, dtₙ₋₂)
end

function alg_cache(
        alg::QNDF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = zero(inv((1 - alg.kappa))), 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    D = Array{typeof(u)}(undef, 1, 2)
    D2 = Array{typeof(u)}(undef, 1, 3)
    R = fill(zero(t), 2, 2)
    U = fill(zero(t), 2, 2)

    D[1] = zero(u)
    D[2] = zero(u)
    D2[1] = zero(u)
    D2[2] = zero(u)
    D2[3] = zero(u)

    U!(2, U)

    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)
    utilde = zero(u)
    uprev2 = zero(u)
    uprev3 = zero(u)
    dtₙ₋₁ = zero(dt)
    dtₙ₋₂ = zero(dt)

    # Migrated `utilde`/`atmp` plus the `nlsolver.tmp` donor alias — zero new
    # arrays vs. the historical cache.
    tmp_cache = TmpCache(nlsolver.tmp, utilde, atmp, nothing, nothing, nothing)

    return QNDF2Cache(
        uprev2, uprev3, fsalfirst, D, D2, R, U, tmp_cache,
        nlsolver, dtₙ₋₁, dtₙ₋₂, alg.step_limiter!
    )
end

@cache mutable struct QNDFConstantCache{
        MO,
        N,
        coefType,
        UType,
        RUType,
        dtType,
        EEstType,
        gammaType,
    } <: OrdinaryDiffEqConstantCache
    nlsolver::N
    U::UType
    R::RUType
    RU::RUType
    D::coefType
    prevD::coefType
    prevorder::Int
    order::Int
    max_order::Val{MO}
    dtprev::dtType
    nconsteps::Int ##Successful Consecutive Step with the same step size
    consfailcnt::Int #Consecutive failed steps count
    EEst1::EEstType #Error Estimator for k-1 order
    EEst2::EEstType #Error Estimator for k+1 order
    γₖ::gammaType
end

function alg_cache(
        alg::QNDF{MO}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits,
    } where {MO}
    max_order = MO
    γ, c = Int64(1) // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    dtprev = one(dt)
    if u isa Number
        D = zeros(uEltypeNoUnits, max_order + 2)
        prevD = zeros(uEltypeNoUnits, max_order + 2)
    else
        D = [zero(u) .* zero(uEltypeNoUnits) for _ in 1:(max_order + 2)]
        prevD = [zero(u) .* zero(uEltypeNoUnits) for _ in 1:(max_order + 2)]
    end
    EEst1 = tTypeNoUnits(1)
    EEst2 = tTypeNoUnits(1)

    U = zeros(tTypeNoUnits, max_order, max_order)
    for r in 1:max_order
        U[1, r] = -r
        for j in 2:max_order
            U[j, r] = U[j - 1, r] * ((j - 1) - r) / j
        end
    end
    R = zeros(tTypeNoUnits, max_order, max_order)
    RU = zeros(tTypeNoUnits, max_order, max_order)

    γₖ = ntuple(k -> sum(tTypeNoUnits(Int64(1) // j) for j in 1:k), Val(max_order))

    return QNDFConstantCache(
        nlsolver, U, R, RU, D, prevD, 1, 1, Val(max_order), dtprev, 0, 0, EEst1,
        EEst2, γₖ
    )
end

@cache mutable struct QNDFCache{
        MO, UType, RUType, rateType, N, coefType, dtType, EEstType,
        gammaType, uType, uNoUnitsType, StepLimiter, TmpC <: TmpCache,
    } <:
    BDFMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    dd::uType
    utildem1::uType
    utildep1::uType
    ϕ::uType
    u₀::uType
    nlsolver::N
    U::UType
    R::RUType
    RU::RUType
    D::coefType
    Dtmp::coefType
    # Unified scratch: the former inline `tmp2` (state scratch, mass-matrix
    # branch) lives in the `tmp` slot, `utilde` in `tmp2`, and `atmp` in
    # `atmp`. The k±1 error buffers (`utildem1`/`utildep1`,
    # `atmpm1`/`atmpp1`) keep their dedicated fields. Rate slots stay opted
    # out (no `preallocate_initdt_buffers` knob on QNDF):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    prevD::coefType
    order::Int
    prevorder::Int
    max_order::Val{MO}
    dtprev::dtType
    nconsteps::Int ##Successful consecutive step with the same step size
    consfailcnt::Int #Consecutive failed steps count
    EEst1::EEstType #Error Estimator for k-1 order
    EEst2::EEstType #Error Estimator for k+1 order
    γₖ::gammaType
    atmpm1::uNoUnitsType
    atmpp1::uNoUnitsType
    dense::Vector{uType}
    step_limiter!::StepLimiter
end

@truncate_stacktrace QNDFCache 1

function alg_cache(
        alg::QNDF{MO}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits,
    } where {MO}
    max_order = MO
    γ, c = Int64(1) // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)
    dd = zero(u)
    utilde = zero(u)
    utildem1 = zero(u)
    utildep1 = zero(u)
    ϕ = zero(u)
    u₀ = zero(u)
    dtprev = one(dt)
    D = [zero(similar(u, uEltypeNoUnits)) for _ in 1:(max_order + 2)]
    Dtmp = [zero(similar(u, uEltypeNoUnits)) for _ in 1:(max_order + 2)]
    prevD = [zero(similar(u, uEltypeNoUnits)) for _ in 1:(max_order + 2)]
    atmp = zero(similar(u, uEltypeNoUnits))
    atmpm1 = zero(similar(u, uEltypeNoUnits))
    atmpp1 = zero(similar(u, uEltypeNoUnits))
    tmp2 = zero(u)
    EEst1 = tTypeNoUnits(1)
    EEst2 = tTypeNoUnits(1)

    U = zeros(tTypeNoUnits, max_order, max_order)
    for r in 1:max_order
        U[1, r] = -r
        for j in 2:max_order
            U[j, r] = U[j - 1, r] * ((j - 1) - r) / j
        end
    end

    R = zeros(tTypeNoUnits, max_order, max_order)
    RU = zeros(tTypeNoUnits, max_order, max_order)
    γₖ = ntuple(k -> sum(tTypeNoUnits(Int64(1) // j) for j in 1:k), Val(max_order))

    dense = [zero(u) for _ in 1:max_order]

    # Migrated fields only (`tmp2` → tmp, `utilde` → tmp2, `atmp` → atmp) —
    # the array count matches the historical cache exactly.
    tmp_cache = TmpCache(tmp2, utilde, atmp, nothing, nothing, nothing)

    return QNDFCache(
        u, uprev, fsalfirst, dd, utildem1, utildep1, ϕ, u₀, nlsolver, U, R, RU, D, Dtmp,
        tmp_cache, prevD, 1, 1, Val(max_order), dtprev, 0, 0, EEst1, EEst2, γₖ,
        atmpm1, atmpp1, dense, alg.step_limiter!
    )
end

@cache mutable struct MEBDF2Cache{uType, rateType, N, TmpC <: TmpCache} <:
    BDFMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    # Unified scratch: the former inline `tmp2` (stage-3 state scratch) lives
    # in the `tmp` slot and `atmp` in `atmp`; no `utilde` existed so `tmp2`
    # stays `nothing`. Rate slots stay opted out (no
    # `preallocate_initdt_buffers` knob on MEBDF2):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    nlsolver::N
end

function alg_cache(
        alg::MEBDF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    tmp2 = zero(u)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    # Migrated fields only (`tmp2` → tmp, `atmp` → atmp) — the array count
    # matches the historical cache exactly.
    tmp_cache = TmpCache(tmp2, nothing, atmp, nothing, nothing, nothing)

    return MEBDF2Cache(u, uprev, uprev2, fsalfirst, z₁, z₂, tmp_cache, nlsolver)
end

mutable struct MEBDF2ConstantCache{N} <: OrdinaryDiffEqConstantCache
    nlsolver::N
end

function alg_cache(
        alg::MEBDF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    return MEBDF2ConstantCache(nlsolver)
end

@cache mutable struct FBDFConstantCache{
        MO, N, tsType, tType, uType, uuType, coeffType,
        EEstType, rType, wType, fdWeightsType, staldType,
    } <:
    OrdinaryDiffEqConstantCache
    nlsolver::N
    ts::tsType
    ts_tmp::tsType
    t_old::tType
    u_history::uuType
    order::Int
    prev_order::Int
    u_corrector::uType
    bdf_coeffs::coeffType
    max_order::Val{MO}
    nconsteps::Int # consecutive success steps
    consfailcnt::Int #consecutive failed step counts
    qwait::Int # countdown to next order change consideration (CVODE-style)
    terkm2::EEstType
    terkm1::EEstType
    terk::EEstType
    terkp1::EEstType
    r::rType
    weights::wType
    iters_from_event::Int
    fd_weights::fdWeightsType
    stald::staldType
end

function alg_cache(
        alg::FBDF{MO}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {
        uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits,
    } where {MO}
    γ, c = Int64(1) // 1, 1
    max_order = MO
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )
    bdf_coeffs = _make_bdf_coeffs_fbdf()
    ts = zero(Vector{typeof(t)}(undef, max_order + 2)) #ts is the successful past points, it will be updated after successful step
    ts_tmp = similar(ts)

    if u isa Number
        u_history = zeros(eltype(u), max_order + 2)
        u_corrector = zeros(eltype(u), max_order + 2)
    else
        u_history = [zero(u) for _ in 1:(max_order + 2)]
        u_corrector = [zero(u) for _ in 1:(max_order + 2)]
    end
    order = 1
    prev_order = 1
    terkm2 = tTypeNoUnits(1)
    terkm1 = tTypeNoUnits(1)
    terk = tTypeNoUnits(1)
    terkp1 = tTypeNoUnits(1)
    r = zero(Vector{typeof(t)}(undef, max_order + 2))
    weights = zero(Vector{typeof(t)}(undef, max_order + 2))
    weights[1] = 1
    nconsteps = 0
    consfailcnt = 0
    qwait = 3 # order + 2, matching nconsteps >= order + 2 for failure-free runs
    t_old = zero(t)
    iters_from_event = 0

    T_stald = real(uBottomEltypeNoUnits)
    stald = StabilityLimitDetectionState(
        T_stald;
        enabled = alg.stald,
        rrcut = alg.stald_rrcut,
        vrrtol = alg.stald_vrrtol,
        vrrt2 = alg.stald_vrrt2,
        sqtol = alg.stald_sqtol,
        rrtol = alg.stald_rrtol,
        tiny = alg.stald_tiny,
    )

    fd_weights = zeros(typeof(t), max_order + 1, max_order + 1)

    return FBDFConstantCache(
        nlsolver, ts, ts_tmp, t_old, u_history, order, prev_order,
        u_corrector, bdf_coeffs, Val(MO), nconsteps, consfailcnt, qwait, terkm2,
        terkm1, terk, terkp1, r, weights, iters_from_event, fd_weights, stald
    )
end

@cache mutable struct FBDFCache{
        MO, N, rateType, tsType, tType, uType, uuType,
        coeffType, EEstType, rType, wType, StepLimiter, fdWeightsType, staldType,
        TmpC <: TmpCache,
    } <:
    BDFMutableCache
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
    nconsteps::Int # consecutive success steps
    consfailcnt::Int #consecutive failed step counts
    qwait::Int # countdown to next order change consideration (CVODE-style)
    # Unified scratch: the former inline `tmp` (corrector RHS, write-first
    # every step) lives in the `tmp` slot and `atmp` in `atmp`; no `utilde`
    # existed so `tmp2` stays `nothing`. The error-estimate buffers
    # `terk_tmp`/`terkp1_tmp` keep their dedicated fields. Rate slots stay
    # opted out (no `preallocate_initdt_buffers` knob on FBDF):
    # `TmpCache{uType, Nothing, uNoUnitsType, Nothing}`.
    tmp_cache::TmpC
    terkm2::EEstType
    terkm1::EEstType
    terk::EEstType #terk corresponds to hᵏyᵏ(tₙ₊₁)
    terkp1::EEstType
    terk_tmp::uType
    terkp1_tmp::uType
    r::rType
    weights::wType #weights of Lagrangian formula
    equi_ts::tsType
    iters_from_event::Int
    dense::Vector{uType}
    step_limiter!::StepLimiter
    fd_weights::fdWeightsType
    stald::staldType
end

@truncate_stacktrace FBDFCache 1

function alg_cache(
        alg::FBDF{MO}, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {
        MO, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits,
    }
    γ, c = Int64(1) // 1, 1
    fsalfirst = zero(rate_prototype)
    max_order = MO
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    bdf_coeffs = _make_bdf_coeffs_fbdf()
    ts = Vector{typeof(t)}(undef, max_order + 2) #ts is the successful past points, it will be updated after successful step
    u_history = [zero(u) for _ in 1:(max_order + 2)]
    order = 1
    prev_order = 1
    u_corrector = [zero(u) for _ in 1:(max_order + 2)]
    recursivefill!(ts, zero(t))
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
    qwait = 3 # order + 2, matching nconsteps >= order + 2 for failure-free runs
    t_old = zero(t)
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, zero(uEltypeNoUnits))
    u₀ = similar(u)
    equi_ts = similar(ts)
    tmp = similar(u)
    ts_tmp = similar(ts)
    iters_from_event = 0

    dense = [zero(u) for _ in 1:(2 * (max_order + 1))]  # first half for integrator.k, second half as scratch

    fd_weights = zeros(typeof(t), max_order + 1, max_order + 1)
    T_stald = real(uBottomEltypeNoUnits)
    stald = StabilityLimitDetectionState(
        T_stald;
        enabled = alg.stald,
        rrcut = alg.stald_rrcut,
        vrrtol = alg.stald_vrrtol,
        vrrt2 = alg.stald_vrrt2,
        sqtol = alg.stald_sqtol,
        rrtol = alg.stald_rrtol,
        tiny = alg.stald_tiny,
    )

    # Migrated fields only (`tmp` → tmp, `atmp` → atmp) — the array count
    # matches the historical cache exactly.
    tmp_cache = TmpCache(tmp, nothing, atmp, nothing, nothing, nothing)

    return FBDFCache(
        fsalfirst, nlsolver, ts, ts_tmp, t_old, u_history, order, prev_order,
        u_corrector, u₀, bdf_coeffs, Val(MO), nconsteps, consfailcnt, qwait, tmp_cache,
        terkm2, terkm1, terk, terkp1, terk_tmp, terkp1_tmp, r, weights, equi_ts,
        iters_from_event, dense, alg.step_limiter!, fd_weights, stald
    )
end
