module OrdinaryDiffEqStabilizedIRK

import OrdinaryDiffEqCore: alg_order, alg_maximum_order,
                           calculate_residuals!,
                           beta2_default, beta1_default, gamma_default, issplit,
                           initialize!, perform_step!, unwrap_alg,
                           calculate_residuals, fac_default_gamma,
                           OrdinaryDiffEqAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqAdaptiveAlgorithm, @SciMLMessage,
                           OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                           alg_cache, _unwrap_val, DEFAULT_PRECS, @cache,
                           _reshape, _vec, full_cache, get_fsalfirstlast,
                           generic_solver_docstring, _bool_to_ADType, _process_AD_choice

using OrdinaryDiffEqDifferentiation: dolinsolve, update_W!
using OrdinaryDiffEqNonlinearSolve: NLNewton, nlsolve!, isnewton, build_nlsolver,
                                    markfirststage!, du_alias_or_new, get_W

using OrdinaryDiffEqStabilizedRK: ESERK4, ESERK5, RKC, SERK2

using FastBroadcast, MuladdMacro, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA
import OrdinaryDiffEqCore
import ADTypes: AutoForwardDiff, AbstractADType

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("irkc_utils.jl")
include("irkc_caches.jl")
include("irkc_perform_step.jl")

export IRKC

end
