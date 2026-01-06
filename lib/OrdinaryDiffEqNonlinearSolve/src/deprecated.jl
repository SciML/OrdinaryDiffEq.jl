# Backwards compatibility layer for build_nlsolver
#
# Downstream packages (StochasticDiffEq, etc.) that CALL build_nlsolver with the
# OLD signature (without verbose) will be forwarded to the NEW signature.
#
# Note: Unlike alg_cache, downstream packages do not DEFINE their own build_nlsolver
# methods - they only CALL it. So we only need the OLD->NEW fallback, not NEW->OLD.

# OLD signature -> NEW signature (for old callers)
# The OLD signature has: (alg, u, uprev, p, t, dt, f, rate_prototype, Types..., γ, c, iip)
# We forward to the NEW signature adding α=1 and verbose=NonlinearVerbosity()
function build_nlsolver(
        alg, u, uprev, p, t, dt, f::F, rate_prototype,
        ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits},
        ::Type{tTypeNoUnits}, γ, c,
        iip
    ) where {F, uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    return build_nlsolver(
        alg, u, uprev, p, t, dt, f, rate_prototype, uEltypeNoUnits,
        uBottomEltypeNoUnits,
        tTypeNoUnits, γ, c, 1, iip, NonlinearVerbosity()
    )
end
