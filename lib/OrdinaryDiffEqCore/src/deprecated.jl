# Backwards compatibility layer for alg_cache
#
# This file provides OLD -> NEW fallbacks for the verbose parameter:
# When old callers don't pass verbose, add ODEVerbosity()
#
# NOTE: NEW -> OLD fallbacks (for calling downstream packages that don't have
# the verbose parameter yet) are intentionally omitted because they cause method
# ambiguity with algorithm-specific methods. Downstream packages must update
# their alg_cache signatures to include the verbose parameter.

# OLD signature -> NEW signature (for old callers using updated OrdinaryDiffEq algorithms)
# These add the default ODEVerbosity() when called without verbose parameter.
function alg_cache(alg, u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{true})
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(true), ODEVerbosity())
end

function alg_cache(alg, u, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits, tTypeNoUnits, uprev, uprev2, f, t,
        dt, reltol, p, calck,
        ::Val{false})
    alg_cache(alg, u, rate_prototype, uEltypeNoUnits, uBottomEltypeNoUnits,
        tTypeNoUnits, uprev, uprev2, f, t, dt, reltol, p, calck, Val(false), ODEVerbosity())
end
