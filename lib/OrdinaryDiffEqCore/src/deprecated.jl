# Backwards compatibility layer for alg_cache
#
# This provides bidirectional fallbacks to support both:
# 1. OLD callers (using old signature without verbose) calling NEW algorithms
# 2. NEW callers (OrdinaryDiffEqCore.__init with verbose) calling OLD algorithms
#
# For downstream packages (DelayDiffEq, StochasticDiffEq, etc.) that haven't updated
# their alg_cache signatures yet, the NEW->OLD fallback allows them to work without changes.
# Once they add the verbose parameter, they can use the full verbosity system.

# OLD signature -> NEW signature (for old callers using updated OrdinaryDiffEq algorithms)
function alg_cache(alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(true), ODEVerbosity())
end

function alg_cache(alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false), ODEVerbosity())
end

# NEW signature -> OLD signature (for downstream packages that haven't added verbose yet)
# These fallbacks strip the verbose argument and call the old signature.
# They only activate when no more specific method (for a particular algorithm type) exists.
function alg_cache(alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true}, verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(true))
end

function alg_cache(alg, u, rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false}, verbose) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false))
end
