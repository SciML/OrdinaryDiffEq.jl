module OrdinaryDiffEqFIRK

import OrdinaryDiffEqCore: unwrap_alg,
    default_controller, PredictiveController, PIController,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    alg_cache, @threaded, isthreaded,
    constvalue,
    differentiation_rk_docstring, trivial_limiter!,
    qmax_default, alg_adaptive_order,
    step_accept_controller!, step_reject_controller!,
    alg_can_repeat_jac, NewtonAlgorithm,
    get_current_adaptive_order, get_fsalfirstlast,
    get_current_alg_order,
    isfirk,
    _fixup_ad, perform_step!,
    LinearAliasSpecifier, set_discontinuity,
    Convergence, FastConvergence, NLStatus,
    VerySlowConvergence, Divergence, get_new_W_γdt_cutoff
import SciMLBase
import SciMLBase: alg_order, _vec, _reshape, _unwrap_val,
    UDerivativeWrapper, UJacobianWrapper, value,
    LinearProblem, get_tmp_cache, init
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
using MuladdMacro: @muladd
using RecursiveArrayTools: recursivefill!
import Polyester
using LinearAlgebra: I, UniformScaling, mul!, lu, dot
import LinearSolve
import FastBroadcast: @..
import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: _ode_interpolant, _ode_interpolant!, _ode_addsteps!,
    has_stiff_interpolation
import FastPower: fastpower
using OrdinaryDiffEqDifferentiation: build_J_W, build_jac_config,
    calc_J!, dolinsolve, calc_J,
    islinearfunction
import ADTypes: AutoForwardDiff
import SciMLOperators
import SciMLOperators: AbstractSciMLOperator
# Load-bearing runtime dependency: provides the nonlinear-solver machinery the FIRK
# integrators dispatch into. The convergence-state names FIRK uses are owned by
# OrdinaryDiffEqCore (imported above), so this is a module-only import to keep the
# dependency from registering as stale.
import OrdinaryDiffEqNonlinearSolve

import OrdinaryDiffEqCore: PredictiveControllerCache

@static if Base.pkgversion(OrdinaryDiffEqCore) >= v"3.10"
    @eval begin
        import OrdinaryDiffEqCore: get_current_qmax
    end
else
    @eval begin
        # Fallback for older OrdinaryDiffEqCore: no first-step qmax behavior
        @inline get_current_qmax(integrator, qmax) = qmax
    end
end

using Reexport: @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("firk_caches.jl")
include("firk_tableaus.jl")
include("firk_perform_step.jl")
include("firk_interpolants.jl")
include("firk_addsteps.jl")
include("integrator_interface.jl")

export RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau, GaussLegendre

end
