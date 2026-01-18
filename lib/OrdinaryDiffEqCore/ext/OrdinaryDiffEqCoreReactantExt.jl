module OrdinaryDiffEqCoreReactantExt

using OrdinaryDiffEqCore
import Reactant
using Reactant: @reactant_overlay, TracedRArray, TracedRNumber

# Reactant cannot convert Rational{Int64} to traced types.
# Override functions that return Rationals to return Float64 instead.

@reactant_overlay function OrdinaryDiffEqCore.qmin_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 0.2 : 0
end

@reactant_overlay function OrdinaryDiffEqCore.beta2_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 2 / (5 * OrdinaryDiffEqCore.alg_order(alg)) : 0
end

@reactant_overlay function OrdinaryDiffEqCore.beta1_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm},
    beta2
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 7 / (10 * OrdinaryDiffEqCore.alg_order(alg)) : 0
end

@reactant_overlay function OrdinaryDiffEqCore.gamma_default(
    alg::Union{OrdinaryDiffEqCore.OrdinaryDiffEqAlgorithm, OrdinaryDiffEqCore.DAEAlgorithm}
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 0.9 : 0
end

@reactant_overlay function OrdinaryDiffEqCore.qsteady_max_default(
    alg::OrdinaryDiffEqCore.OrdinaryDiffEqAdaptiveImplicitAlgorithm
)
    return 1.2
end

@reactant_overlay function OrdinaryDiffEqCore.qsteady_max_default(
    alg::OrdinaryDiffEqCore.OrdinaryDiffEqImplicitAlgorithm
)
    return OrdinaryDiffEqCore.isadaptive(alg) ? 1.0 : 0
end

# NOTE: Full Reactant support for ODE solving also requires Reactant-aware
# handling in DiffEqBase and RecursiveArrayTools to avoid undefined array
# references during solution saving. The functions copyat_or_push! and
# recursivecopy use `similar()` to pre-allocate arrays with undefined elements,
# which Reactant cannot trace. This would need upstream changes to:
# - RecursiveArrayTools: Add Reactant extension with push-based alternatives
# - DiffEqBase: Add Reactant-aware saving path in resize_non_user_cache!

end
