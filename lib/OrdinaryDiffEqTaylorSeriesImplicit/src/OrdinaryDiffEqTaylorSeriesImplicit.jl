module OrdinaryDiffEqTaylorSeriesImplicit

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order, alg_autodiff, calculate_residuals!,
    initialize!, perform_step!, unwrap_alg,
    calculate_residuals, alg_extrapolates,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqAdaptiveAlgorithm,
    alg_cache, _vec, _reshape, @cache,
    constvalue, _unwrap_val,
    trivial_limiter!, get_fsalfirstlast, generic_solver_docstring,
    Convergence, NLStatus, VerySlowConvergence, Divergence
import OrdinaryDiffEqCore
using FastBroadcast: @..
using RecursiveArrayTools: recursivefill!
using SciMLBase: SciMLBase, LinearProblem, LinearAliasSpecifier, init
using OrdinaryDiffEqDifferentiation: dolinsolve
import ADTypes: AutoForwardDiff
import FunctionWrappers: FunctionWrapper
import LinearSolve
import ArrayInterface

include("algorithms.jl")
include("alg_utils.jl")
include("TaylorSeries_caches.jl")
include("TaylorSeries_perform_step.jl")

export ImplicitTaylor
end
