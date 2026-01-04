@cache struct PDIRK44Cache{uType, rateType, N, TabType} <: OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    k1::Array{rateType}
    k2::Array{rateType}
    nlsolver::N
    tab::TabType
end

# Non-FSAL
get_fsalfirstlast(cache::PDIRK44Cache, u) = (nothing, nothing)

struct PDIRK44ConstantCache{N, TabType} <: OrdinaryDiffEqConstantCache
    nlsolver::N
    tab::TabType
end

struct PDIRK44Tableau{T, T2}
    γs::SVector{2, T2}
    cs::SVector{4, T2}
    α1::SVector{2, T}
    α2::SVector{2, T}
    b1::T
    b2::T
    b3::T
    b4::T
end

function PDIRK44Tableau(T, T2)
    γ1 = convert(T2, 1 // 2)
    γ2 = convert(T2, 2 // 3)
    γs = SVector(γ1, γ2)
    c1 = convert(T2, 1 // 2)
    c2 = convert(T2, 2 // 3)
    c3 = convert(T2, 1 // 2)
    c4 = convert(T2, 1 // 3)
    cs = SVector(c1, c2, c3, c4)
    α11 = convert(T, -5 // 2)
    α12 = convert(T, -5 // 3)
    α1 = SVector(α11, α12)
    α21 = convert(T, 5 // 2)
    α22 = convert(T, 4 // 3)
    α2 = SVector(α21, α22)
    b1 = convert(T, -1 // 1)
    b2 = convert(T, -1 // 1)
    b3 = convert(T, 3 // 2)
    b4 = convert(T, 3 // 2)
    return PDIRK44Tableau{T, T2}(γs, cs, α1, α2, b1, b2, b3, b4)
end

function alg_cache(
        alg::PDIRK44, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1.0, 1.0
    if alg.threading
        nlsolver1 = build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype,
            uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c,
            Val(true)
        )
        nlsolver2 = build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype,
            uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c,
            Val(true)
        )
        nlsolver = [nlsolver1, nlsolver2]
    else
        _nlsolver = build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype,
            uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c,
            Val(true)
        )
        nlsolver = [_nlsolver]
    end
    tab = PDIRK44Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    k1 = [zero(rate_prototype) for i in 1:2]
    k2 = [zero(rate_prototype) for i in 1:2]
    return PDIRK44Cache(u, uprev, k1, k2, nlsolver, tab)
end

function alg_cache(
        alg::PDIRK44, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    γ, c = 1.0, 1.0
    if alg.threading
        nlsolver1 = build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype,
            uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c,
            Val(false)
        )
        nlsolver2 = build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype,
            uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c,
            Val(false)
        )
        nlsolver = [nlsolver1, nlsolver2]
    else
        _nlsolver = build_nlsolver(
            alg, u, uprev, p, t, dt, f, rate_prototype,
            uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits, γ, c,
            Val(false)
        )
        nlsolver = [_nlsolver]
    end
    tab = PDIRK44Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    return PDIRK44ConstantCache(nlsolver, tab)
end
