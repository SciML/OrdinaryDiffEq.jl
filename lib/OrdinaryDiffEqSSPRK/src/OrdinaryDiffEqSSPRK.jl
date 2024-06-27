module OrdinaryDiffEqSSPRK

import OrdinaryDiffEq: alg_order, alg_adaptive_order, calculate_residuals!,
                       beta2_default, beta1_default, gamma_default,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, ssp_coefficient,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new, ArrayFuse
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

include("alg_utils.jl")
include("algorithms.jl")
include("ssprk_caches.jl")
include("ssprk_perform_step.jl")

export SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
       SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H, 
       SSPRK33, SHLDDRK_2N, KYKSSPRK42, SHLDDRK52

end