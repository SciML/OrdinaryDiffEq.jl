module OrdinaryDiffEqExtrapolation

import OrdinaryDiffEqCore: alg_order, alg_maximum_order, get_current_adaptive_order,
    get_current_alg_order, calculate_residuals!,
    accept_step_controller,
    beta2_default, beta1_default, gamma_default,
    initialize!, perform_step!, @cache, unwrap_alg,
    isthreaded, isadaptive, PIController,
    step_accept_controller!, calculate_residuals,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    reset_alg_dependent_opts!, AbstractController,
    step_accept_controller!, step_reject_controller!,
    OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqAdaptiveImplicitAlgorithm,
    alg_cache, CompiledFloats, @threaded, stepsize_controller!,
    DEFAULT_PRECS, full_cache, qmin_default, qmax_default,
    constvalue, PolyesterThreads, Sequential, BaseThreads,
    _digest_beta1_beta2, timedepentdtmin, _unwrap_val,
    _reshape, _vec, get_fsalfirstlast, generic_solver_docstring,
    differentiation_rk_docstring, _bool_to_ADType,
    _process_AD_choice, LinearAliasSpecifier, @SciMLMessage, Minimal

using FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools, LinearSolve
import OrdinaryDiffEqCore
import FastPower
import OrdinaryDiffEqDifferentiation: TimeDerivativeWrapper, UDerivativeWrapper, calc_J,
    WOperator, TimeGradientWrapper, UJacobianWrapper,
    build_grad_config,
    build_jac_config, calc_J!, jacobian2W!, dolinsolve
import ADTypes: AutoForwardDiff, AbstractADType

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.3"
    @eval begin
        import OrdinaryDiffEqCore: AbstractLegacyController, default_controller_v7,
            legacy_default_controller, setup_controller_cache, AbstractControllerCache
    end
end

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("extrapolation_caches.jl")
include("extrapolation_perform_step.jl")

@inline function SciMLBase.get_tmp_cache(
        integrator,
        alg::OrdinaryDiffEqImplicitExtrapolationAlgorithm,
        cache::OrdinaryDiffEqMutableCache
    )
    return (cache.tmp, cache.utilde)
end

export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
    ImplicitEulerExtrapolation,
    ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
    ImplicitEulerBarycentricExtrapolation

end
