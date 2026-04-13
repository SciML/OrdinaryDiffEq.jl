#Backwards compat for StochasticDiffEq
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
