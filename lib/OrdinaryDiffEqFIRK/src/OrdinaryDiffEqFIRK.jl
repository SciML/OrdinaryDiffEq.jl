module OrdinaryDiffEqFIRK

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals,
                       OrdinaryDiffEqAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new,
                       explicit_rk_docstring, trivial_limiter!,
                       _ode_interpolant!, _ode_addsteps!, AbstractController,
                       NLStatus, qmax_default, alg_adaptive_order, DEFAULT_PRECS,
                       UJacobianWrapper, build_J_W, build_jac_config, UDerivativeWrapper,
                       Convergence, calc_J!, dolinsolve, FastConvergence, calc_J,
                       stepsize_controller!, islinearfunction, step_accept_controller!, step_reject_controller!
using MuladdMacro, DiffEqBase, RecursiveArrayTools
using SciMLOperators: AbstractSciMLOperator
using LinearAlgebra: I, UniformScaling, mul!, lu
import LinearSolve
import FastBroadcast: @..
include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("firk_caches.jl")
include("firk_tableaus.jl")
include("firk_perform_step.jl")
include("integrator_interface.jl")

export RadauIIA3, RadauIIA5, RadauIIA7

end
