module OrdinaryDiffEqExtrapolation

import OrdinaryDiffEqCore: TmpCache, build_tmp_cache,
    alg_maximum_order, get_current_adaptive_order,
    get_current_alg_order,
    accept_step_controller,
    beta2_default, beta1_default, gamma_default,
    perform_step!, @cache, unwrap_alg,
    isthreaded, PIController,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    reset_alg_dependent_opts!, AbstractController,
    step_accept_controller!, step_reject_controller!,
    OrdinaryDiffEqAdaptiveAlgorithm,
    OrdinaryDiffEqAdaptiveImplicitAlgorithm,
    alg_cache, CompiledFloats, @threaded, stepsize_controller!,
    full_cache, qmin_default,
    constvalue, PolyesterThreads,
    _fixup_ad,
    get_fsalfirstlast, generic_solver_docstring,
    differentiation_rk_docstring,
    LinearAliasSpecifier
import OrdinaryDiffEqCore: default_controller, AbstractControllerCache, setup_controller_cache,
    get_qmin, get_qmax
import OrdinaryDiffEqCore

# Owned by SciMLBase / DiffEqBase / SciMLOperators, re-exported through
# OrdinaryDiffEqCore / OrdinaryDiffEqDifferentiation — import from the owners directly so
# ExplicitImports' owner check is satisfied.
import SciMLBase: alg_order, _unwrap_val, _reshape, _vec,
    TimeDerivativeWrapper, UDerivativeWrapper,
    TimeGradientWrapper, UJacobianWrapper
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!, timedepentdtmin

using FastBroadcast: FastBroadcast, @..
using MuladdMacro: MuladdMacro, @muladd
using RecursiveArrayTools: RecursiveArrayTools, recursivefill!
using LinearSolve: LinearSolve, RFLUFactorization
import FastPower
using SciMLOperators: SciMLOperators, WOperator
import SciMLLogging: @SciMLMessage
import OrdinaryDiffEqDifferentiation: calc_J,
    build_grad_config,
    build_jac_config, calc_J!, jacobian2W!, dolinsolve
import ADTypes: AutoForwardDiff

using Reexport: Reexport, @reexport
@reexport using SciMLBase
using SciMLBase: SciMLBase, LinearProblem, init

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
    return (cache.tmp_cache.tmp, cache.tmp_cache.tmp2)
end

export AitkenNeville, ExtrapolationMidpointDeuflhard, ExtrapolationMidpointHairerWanner,
    ImplicitEulerExtrapolation,
    ImplicitDeuflhardExtrapolation, ImplicitHairerWannerExtrapolation,
    ImplicitEulerBarycentricExtrapolation

end
