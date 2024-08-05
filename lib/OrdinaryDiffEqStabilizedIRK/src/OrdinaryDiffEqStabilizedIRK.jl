module OrdinaryDiffEqStabilizedIRK

import OrdinaryDiffEq: alg_order, alg_maximum_order,
                       calculate_residuals!,
                       beta2_default, beta1_default, gamma_default, issplit,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, fac_default_gamma,
                       OrdinaryDiffEqAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm,
                       OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                       alg_cache, _unwrap_val, DEFAULT_PRECS, @cache,
                       _reshape, _vec

using OrdinaryDiffEq.OrdinaryDiffEqDifferentiation: dolinsolve, update_W!
using OrdinaryDiffEq.OrdinaryDiffEqNonlinearSolve: NLNewton, nlsolve!, isnewton, build_nlsolver, markfirststage!, du_alias_or_new, get_W
using DiffEqBase, FastBroadcast, MuladdMacro, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

include("algorithms.jl")
include("alg_utils.jl")
include("irkc_utils.jl")
include("irkc_caches.jl")
include("irkc_perform_step.jl")

export IRKC

end
