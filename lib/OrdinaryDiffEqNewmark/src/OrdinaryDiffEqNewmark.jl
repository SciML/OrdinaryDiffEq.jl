module OrdinaryDiffEqNewmark

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, @unpack, unwrap_alg,
                           calculate_residuals, alg_extrapolates,
                           OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqNewtonAlgorithm,
                           DEFAULT_PRECS,
                           OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                           constvalue, _unwrap_val, _ode_interpolant,
                           trivial_limiter!, _ode_interpolant!,
                           isesdirk, issplit,
                           ssp_coefficient, get_fsalfirstlast, generic_solver_docstring
using TruncatedStacktraces, MuladdMacro, MacroTools, FastBroadcast, RecursiveArrayTools
using SciMLBase: DynamicalODEFunction
using LinearAlgebra: mul!, I
import OrdinaryDiffEqCore

using OrdinaryDiffEqDifferentiation: UJacobianWrapper, dolinsolve
using OrdinaryDiffEqNonlinearSolve: du_alias_or_new, markfirststage!, build_nlsolver,
                                    nlsolve!, nlsolvefail, isnewton, get_W, set_new_W!,
                                    NLNewton, COEFFICIENT_MULTISTEP

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("newmark_caches.jl")
include("newmark_nlsolve.jl")
include("newmark_perform_step.jl")

export NewmarkBeta

end