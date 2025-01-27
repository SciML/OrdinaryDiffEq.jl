module OrdinaryDiffEqFIRK

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, @unpack, unwrap_alg,
                           calculate_residuals,
                           OrdinaryDiffEqAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, @threaded, isthreaded, PolyesterThreads, 
                           isfsal, full_cache, constvalue, _unwrap_val,
                           differentiation_rk_docstring, trivial_limiter!,
                           _ode_interpolant!, _ode_addsteps!, AbstractController,
                           qmax_default, alg_adaptive_order, DEFAULT_PRECS,
                           stepsize_controller!, step_accept_controller!,
                           step_reject_controller!,
                           PredictiveController, alg_can_repeat_jac, NewtonAlgorithm,
                           fac_default_gamma,
                           get_current_adaptive_order, get_fsalfirstlast,
                           isfirk, generic_solver_docstring, _bool_to_ADType,
                           _process_AD_choice
using MuladdMacro, DiffEqBase, RecursiveArrayTools, Polyester
                           isfirk, generic_solver_docstring
using SciMLOperators: AbstractSciMLOperator
using LinearAlgebra: I, UniformScaling, mul!, lu
import LinearSolve
import FastBroadcast: @..
import OrdinaryDiffEqCore
import FastPower

using OrdinaryDiffEqDifferentiation: UJacobianWrapper, build_J_W, build_jac_config,
                                     UDerivativeWrapper, calc_J!, dolinsolve, calc_J,
                                     islinearfunction
using OrdinaryDiffEqNonlinearSolve: du_alias_or_new, Convergence, FastConvergence, NLStatus,
                                    VerySlowConvergence,
                                    Divergence, get_new_W_γdt_cutoff
import ADTypes: AutoForwardDiff, AbstractADType

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("firk_caches.jl")
include("firk_tableaus.jl")
include("firk_perform_step.jl")
include("integrator_interface.jl")

export RadauIIA3, RadauIIA5, RadauIIA9, AdaptiveRadau

end
