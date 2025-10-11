module OrdinaryDiffEqNewmark

import OrdinaryDiffEqCore: initialize!, perform_step!, @unpack, unwrap_alg,
                           alg_extrapolates, isadaptive, alg_order, 
                           OrdinaryDiffEqAlgorithm, OrdinaryDiffEqNewtonAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqImplicitSecondOrderAlgorithm,
                           OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm,
                           OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                           constvalue, _unwrap_val, _ode_interpolant,
                           trivial_limiter!, _ode_interpolant!,
                           get_fsalfirstlast, generic_solver_docstring,
                           OrdinaryDiffEqCore
using TruncatedStacktraces, MuladdMacro, MacroTools, FastBroadcast, RecursiveArrayTools
using SciMLBase: DynamicalODEFunction
using LinearAlgebra: mul!, I

using NonlinearSolveFirstOrder

using Reexport
@reexport using DiffEqBase

using ConcreteStructs

include("algorithms.jl")
include("alg_utils.jl")
include("newmark_caches.jl")
include("newmark_nlsolve.jl")
include("newmark_perform_step.jl")

export NewmarkBeta

end
