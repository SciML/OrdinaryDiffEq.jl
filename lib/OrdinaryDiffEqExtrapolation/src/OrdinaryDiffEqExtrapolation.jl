module OrdinaryDiffEqExtrapolation

import OrdinaryDiffEqCore: alg_maximum_order, get_current_adaptive_order,
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
    alg_cache, CompiledFloats, stepsize_controller!,
    qmin_default,
    constvalue, Sequential, BaseThreads, PolyesterThreads,
    _fixup_ad,
    get_fsalfirstlast
import OrdinaryDiffEqCore: default_controller, AbstractControllerCache, setup_controller_cache,
    get_qmin, get_qmax
import OrdinaryDiffEqCore

# Owned by SciMLBase / DiffEqBase, re-exported through OrdinaryDiffEqCore /
# OrdinaryDiffEqDifferentiation.
import SciMLBase: alg_order, LinearAliasSpecifier,
    TimeDerivativeWrapper, UDerivativeWrapper,
    TimeGradientWrapper, UJacobianWrapper
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!, timedepentdtmin

using FastBroadcast: FastBroadcast, @..
using MuladdMacro: MuladdMacro, @muladd
using RecursiveArrayTools: RecursiveArrayTools, recursivefill!
using LinearSolve: LinearSolve, GenericLUFactorization
using FastPower: fastpower
using SciMLOperators: SciMLOperators, WOperator
import SciMLLogging: @SciMLMessage
import OrdinaryDiffEqDifferentiation: calc_J,
    build_grad_config,
    build_jac_config, calc_J!, jacobian2W!, dolinsolve
import ADTypes: AutoForwardDiff

using CommonSolve: init
using SciMLBase: SciMLBase, LinearProblem

using Reexport: Reexport, @reexport
# Kept in sync with docs/src/api/reexports.md by test/qa/reexport_tests.jl.
@reexport using SciMLBase: DAEProblem, DiscreteProblem, DynamicalODEProblem,
    EnsembleProblem, ODEProblem, SecondOrderODEProblem, SplitODEProblem, DAEFunction,
    DiscreteFunction, DynamicalODEFunction, ODEFunction, SplitFunction, solve, solve!, init,
    step!, remake, ReturnCode, successful_retcode, DEStats, NLStats, NullParameters,
    AutoSpecialize, FullSpecialize, NoSpecialize, FunctionWrapperSpecialize, CheckInit,
    NoInit, OverrideInit, CallbackSet, ContinuousCallback, DiscreteCallback,
    VectorContinuousCallback, LeftRootFind, RightRootFind, NoRootFind, add_saveat!,
    add_tstop!, auto_dt_reset!, change_t_via_interpolation!, check_error, check_error!,
    first_tstop, get_dt, get_du, get_du!, get_proposed_dt, get_tmp_cache, pop_tstop!,
    reinit!, savevalues!, set_abstol!, set_proposed_dt!, set_reltol!, set_t!, set_u!,
    set_ut!, terminate!, u_modified!, EnsembleAnalysis, EnsembleDistributed, EnsembleSerial,
    EnsembleSplitThreads, EnsembleSummary, EnsembleThreads

include("utils.jl")
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

@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :Sequential, :BaseThreads, :PolyesterThreads))
end

end
