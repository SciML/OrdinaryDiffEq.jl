module OrdinaryDiffEqExtrapolation

import OrdinaryDiffEq: alg_order, alg_maximum_order, get_current_adaptive_order,
                       get_current_alg_order, calculate_residuals!, accept_step_controller,
                       default_controller, beta2_default, beta1_default, gamma_default,
                       initialize!, perform_step!, @unpack, @cache, unwrap_alg, isthreaded,
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


include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("extrapolation_caches.jl")
include("extrapolation_perform_step.jl")

@inline function DiffEqBase.get_tmp_cache(integrator,
        alg::OrdinaryDiffEqImplicitExtrapolationAlgorithm,
        cache::OrdinaryDiffEqMutableCache)
    (cache.tmp, cache.utilde)
end

export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
       ImplicitEulerExtrapolation,
       ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
       ImplicitEulerBarycentricExtrapolation

end
