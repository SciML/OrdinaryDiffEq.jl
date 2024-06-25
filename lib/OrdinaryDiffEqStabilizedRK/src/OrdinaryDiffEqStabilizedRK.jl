module OrdinaryDiffEqStabilizedRK

import OrdinaryDiffEq: alg_order, alg_adaptive_order, calculate_residuals!,
                       beta2_default, beta1_default, gamma_default,
                       fac_default_gamma, has_dtnew_modification,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals,
                       OrdinaryDiffEqAlgorithm, ispredictive,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, calc_dt_propose!,
                       alg_cache, _vec, _reshape, @cache,
                       constvalue, _unwrap_val, du_alias_or_new
using DiffEqBase, FastBroadcast, MuladdMacro, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

include("algorithms.jl")
include("alg_utils.jl")
include("rkc_utils.jl")
include("rkc_caches.jl")
include("rkc_perform_step.jl")
include("rkc_tableaus_serk2.jl")
include("rkc_tableaus_rock4.jl")
include("rkc_tableaus_rock2.jl")
include("rkc_tableaus_eserk5.jl")
include("rkc_tableaus_eserk4.jl")

export ROCK2, ROCK4, RKC, ESERK4, ESERK5, SERK2

end
