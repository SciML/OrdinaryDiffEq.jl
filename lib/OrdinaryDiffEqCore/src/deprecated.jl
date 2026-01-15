# Backwards compat for old callers (without verbose) -> new methods (with verbose)
# e.g., DelayDiffEq calling new OrdinaryDiffEq
function alg_cache(
        alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return alg_cache(
        alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(true), ODEVerbosity()
    )
end

function alg_cache(
        alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return alg_cache(
        alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false), ODEVerbosity()
    )
end

# Backwards compat for new callers (with verbose) -> old methods (without verbose)
# e.g., new OrdinaryDiffEqCore calling old sublibrary versions from registry
# Julia dispatch prefers more specific methods (typed alg), so these generic
# fallbacks only match when the specific method doesn't accept verbose.
function alg_cache(
        alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        v::Val{true}, ::ODEVerbosity
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return alg_cache(
        alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, v
    )
end

function alg_cache(
        alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        v::Val{false}, ::ODEVerbosity
    ) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return alg_cache(
        alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, v
    )
end
