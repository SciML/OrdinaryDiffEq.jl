@cache struct GLEECache{uType, yType, rateType, yRateType, uNoUnitsType, TabType} <:
    OrdinaryDiffEqCore.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tmp::uType
    fsalfirst::rateType
    fsallast::rateType
    ks::Vector{yRateType}
    ytmp::yType
    εloc::yType
    atmp::uNoUnitsType
    tab::TabType
end

OrdinaryDiffEqCore.get_fsalfirstlast(cache::GLEECache, u) =
    (cache.fsalfirst, cache.fsallast)

struct GLEEConstantCache{TabType} <: OrdinaryDiffEqCore.OrdinaryDiffEqConstantCache
    tab::TabType
end

function OrdinaryDiffEqCore.alg_cache(
        alg::AbstractGLEEAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = _glee_tableau_for(
        alg, OrdinaryDiffEqCore.constvalue(uBottomEltypeNoUnits),
        OrdinaryDiffEqCore.constvalue(tTypeNoUnits)
    )
    y = u.x[1]
    yrate = rate_prototype.x[1]
    ks = [zero(yrate) for _ in 1:nstages(tab)]
    return GLEECache(
        u, uprev, zero(u), zero(rate_prototype), zero(rate_prototype), ks,
        zero(y), zero(y), fill!(similar(y, uEltypeNoUnits), 0), tab
    )
end

function OrdinaryDiffEqCore.alg_cache(
        alg::AbstractGLEEAlgorithm, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    tab = _glee_tableau_for(
        alg, OrdinaryDiffEqCore.constvalue(uBottomEltypeNoUnits),
        OrdinaryDiffEqCore.constvalue(tTypeNoUnits)
    )
    return GLEEConstantCache(tab)
end
