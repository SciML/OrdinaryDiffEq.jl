@cache mutable struct CFNLIRK3ConstantCache{N, Tab} <: SDIRKConstantCache
    nlsolver::N
    tab::Tab
end

function alg_cache(
        alg::CFNLIRK3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
        uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = CFNLIRK3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(false), verbose
    )

    return CFNLIRK3ConstantCache(nlsolver, tab)
end

@cache mutable struct CFNLIRK3Cache{uType, rateType, uNoUnitsType, N, Tab, kType} <:
    SDIRKMutableCache
    u::uType
    uprev::uType
    fsalfirst::rateType
    z₁::uType
    z₂::uType
    z₃::uType
    z₄::uType
    k1::kType
    k2::kType
    k3::kType
    k4::kType
    atmp::uNoUnitsType
    nlsolver::N
    tab::Tab
end

function alg_cache(
        alg::CFNLIRK3, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, uprev, uprev2, f, t, dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = CFNLIRK3Tableau(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
    γ, c = tab.γ, tab.c3
    nlsolver = build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, γ, c, Val(true), verbose
    )
    fsalfirst = zero(rate_prototype)

    k1 = zero(u)
    k2 = zero(u)
    k3 = zero(u)
    k4 = zero(u)

    z₁ = zero(u)
    z₂ = zero(u)
    z₃ = zero(u)
    z₄ = nlsolver.z
    atmp = similar(u, uEltypeNoUnits)
    recursivefill!(atmp, false)

    return CFNLIRK3Cache(u, uprev, fsalfirst, z₁, z₂, z₃, z₄, k1, k2, k3, k4, atmp, nlsolver, tab)
end
