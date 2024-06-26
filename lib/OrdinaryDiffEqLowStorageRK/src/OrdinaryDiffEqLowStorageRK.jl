module OrdinaryDiffEqLowStorageRK

import OrdinaryDiffEq: alg_order, alg_adaptive_order, calculate_residuals!,
                       beta2_default, beta1_default, gamma_default,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, ssp_coefficient,
                       OrdinaryDiffEqAlgorithm, ispredictive,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, uses_uprev,
                       default_controller,
                       alg_cache, _vec, _reshape, @cache, isfsal,
                       constvalue, _unwrap_val, du_alias_or_new, ArrayFuse
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA


include("algorithms.jl")
include("alg_utils.jl")
include("low_storage_rk_caches.jl")
include("low_storage_rk_perform_step.jl")

export ORK256, CarpenterKennedy2N54, SHLDDRK64, HSLDDRK64, DGLDDRK73_C, DGLDDRK84_C,
       DGLDDRK84_F, NDBLSRK124, NDBLSRK134, NDBLSRK144, 
       CFRLDDRK64, TSLDDRK74, CKLLSRK43_2, CKLLSRK54_3C, 
       CKLLSRK95_4S, CKLLSRK95_4C, CKLLSRK95_4M,
       CKLLSRK54_3C_3R, CKLLSRK54_3M_3R, CKLLSRK54_3N_3R, CKLLSRK85_4C_3R, CKLLSRK85_4M_3R, CKLLSRK85_4P_3R,
       CKLLSRK54_3N_4R, CKLLSRK54_3M_4R, CKLLSRK65_4M_4R, CKLLSRK85_4FM_4R, CKLLSRK75_4M_5R, 
       ParsaniKetchesonDeconinck3S32, ParsaniKetchesonDeconinck3S82,
       ParsaniKetchesonDeconinck3S53, ParsaniKetchesonDeconinck3S173,
       ParsaniKetchesonDeconinck3S94, ParsaniKetchesonDeconinck3S184,
       ParsaniKetchesonDeconinck3S105, ParsaniKetchesonDeconinck3S205,
       RDPK3Sp35, RDPK3SpFSAL35, RDPK3Sp49, RDPK3SpFSAL49, RDPK3Sp510, RDPK3SpFSAL510,
       KYK2014DGSSPRK_3S2, RK46NL
end
