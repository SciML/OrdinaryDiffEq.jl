# IMEX Multistep methods
abstract type IMEXMutableCache <: OrdinaryDiffEqMutableCache end
function get_fsalfirstlast(cache::IMEXMutableCache, u)
    return (cache.fsalfirst, du_alias_or_new(cache.nlsolver, cache.fsalfirst))
end

# CNAB2

@cache mutable struct CNAB2ConstantCache{rateType, N, uType, tType} <:
    OrdinaryDiffEqConstantCache
    k2::rateType
    nlsolver::N
    uprev3::uType
    tprev2::tType
end

@cache mutable struct CNAB2Cache{
        uType, rateType, N, tType, TmpC <: TmpCache,
    } <: IMEXMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k1::rateType
    k2::rateType
    du₁::rateType
    nlsolver::N
    uprev3::uType
    tprev2::tType
    # Unified scratch surface. This cache has no inline scratch of its own
    # (state scratch lives in the nlsolver, which is Newton state and off
    # limits), so the state/unit-less slots are opted out (`nothing`). The rate
    # slots donor-alias `k1`/`du₁`, which are dead between steps — both are
    # fully rewritten at the top of every `perform_step!` before being read and
    # never feed dense output — so `initdt` runs allocation-free on the rate
    # side with zero new arrays.
    tmp_cache::TmpC
end

function alg_cache(
        alg::CNAB2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 2, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    k2 = rate_prototype
    uprev3 = u
    tprev2 = t

    return CNAB2ConstantCache(k2, nlsolver, uprev3, tprev2)
end

function alg_cache(
        alg::CNAB2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 2, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    du₁ = zero(rate_prototype)
    uprev3 = zero(u)
    tprev2 = t

    # Footprint-neutral TmpCache: no slot allocates. `k1`/`du₁` are recomputed
    # at the top of every step before any read (dead between steps, never read
    # by interpolation), so they are legal donors for the rate slots. The
    # state/unit-less slots never existed inline here and stay `nothing`.
    tmp_cache = TmpCache(nothing, nothing, nothing, nothing, k1, du₁)

    return CNAB2Cache(
        u, uprev, uprev2, fsalfirst, k1, k2, du₁, nlsolver, uprev3, tprev2, tmp_cache
    )
end

# CNLF2

@cache mutable struct CNLF2ConstantCache{rateType, N, uType, tType} <:
    OrdinaryDiffEqConstantCache
    k2::rateType
    nlsolver::N
    uprev2::uType
    uprev3::uType
    tprev2::tType
end

@cache mutable struct CNLF2Cache{
        uType, rateType, N, tType, TmpC <: TmpCache,
    } <: IMEXMutableCache
    u::uType
    uprev::uType
    uprev2::uType
    fsalfirst::rateType
    k1::rateType
    k2::rateType
    du₁::rateType
    nlsolver::N
    uprev3::uType
    tprev2::tType
    # Unified scratch surface; same layout as `CNAB2Cache`: state/unit-less
    # slots opted out, rate slots donor-aliased to `k1`/`du₁`. `k1` is never
    # read by CNLF2's step at all and `du₁` is rewritten at the top of every
    # `perform_step!` before being read; neither feeds dense output. (`k2` and
    # `uprev2` are multistep history — read before written — so they are NOT
    # donors.)
    tmp_cache::TmpC
end

function alg_cache(
        alg::CNLF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    k2 = rate_prototype
    uprev2 = u
    uprev3 = u
    tprev2 = t

    return CNLF2ConstantCache(k2, nlsolver, uprev2, uprev3, tprev2)
end

function alg_cache(
        alg::CNLF2, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1 // 1, 1
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    k1 = zero(rate_prototype)
    k2 = zero(rate_prototype)
    du₁ = zero(rate_prototype)
    uprev2 = zero(u)
    uprev3 = zero(u)
    tprev2 = t

    # Footprint-neutral TmpCache: no slot allocates. `k1` is never read by this
    # method's step and `du₁` is recomputed each step before any read, so both
    # are legal rate-slot donors (see the struct comment).
    tmp_cache = TmpCache(nothing, nothing, nothing, nothing, k1, du₁)

    return CNLF2Cache(
        u, uprev, uprev2, fsalfirst, k1, k2, du₁, nlsolver, uprev3, tprev2, tmp_cache
    )
end
