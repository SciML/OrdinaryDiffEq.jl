module OrdinaryDiffEqLowStorageRK

import OrdinaryDiffEq: alg_order, alg_maximum_order, get_current_adaptive_order,
                       get_current_alg_order, calculate_residuals!, accept_step_controller,
                       default_controller, beta2_default, beta1_default, gamma_default,
                       initialize!, perform_step!, @unpack, unwrap_alg, isthreaded,
                       step_accept_controller!, calculate_residuals,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       reset_alg_dependent_opts!, AbstractController,
                       step_accept_controller!, step_reject_controller!,
                       OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                       alg_cache, CompiledFloats, @threaded, stepsize_controller!, DEFAULT_PRECS,
                       constvalue, PolyesterThreads, Sequential, BaseThreads,
                       _digest_beta1_beta2, timedepentdtmin, _unwrap_val,
                       TimeDerivativeWrapper, UDerivativeWrapper, calc_J, _reshape, _vec,
                       WOperator, TimeGradientWrapper, UJacobianWrapper, build_grad_config,
                       build_jac_config, calc_J!, jacobian2W!, dolinsolve         
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools, LinearSolve
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

macro cache(expr)
    name = expr.args[2].args[1].args[1]
    fields = [x for x in expr.args[3].args if typeof(x) != LineNumberNode]
    cache_vars = Expr[]
    jac_vars = Pair{Symbol, Expr}[]
    for x in fields
        if x.args[2] == :uType || x.args[2] == :rateType ||
           x.args[2] == :kType || x.args[2] == :uNoUnitsType
            push!(cache_vars, :(c.$(x.args[1])))
        elseif x.args[2] == :DiffCacheType
            push!(cache_vars, :(c.$(x.args[1]).du))
            push!(cache_vars, :(c.$(x.args[1]).dual_du))
        end
    end
    quote
        $(esc(expr))
        $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
    end
end

include("low_storage_rk_caches.jl")
include("low_storage_rk_perform_step.jl")
include("algorithms.jl")
include("alg_utils.jl")

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
