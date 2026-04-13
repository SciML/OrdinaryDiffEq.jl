module StochasticDiffEqImplicit

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!, issplit

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step, supports_regular_jumps, isadaptive,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqNewtonAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    StochasticDiffEqCompositeAlgorithm,
    unwrap_alg,
    @cache

import DiffEqBase: is_diagonal_noise, @.., SplitSDEFunction
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache, has_Wfact

import MuladdMacro: @muladd
import SciMLBase

using OrdinaryDiffEqCore: _vec, _reshape, current_extrapolant, current_extrapolant!,
    isnewton, set_new_W!, get_W
using OrdinaryDiffEqNonlinearSolve: NLSolver, nlsolvefail, nlsolve!, build_nlsolver,
    markfirststage!, NLNewton
import OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqDifferentiation: calc_J, calc_J!, dolinsolve, get_W

using JumpProcesses: JumpProblem
using LinearAlgebra
using StaticArrays
using RecursiveArrayTools
using ForwardDiff, FiniteDiff
import ADTypes

include("tableaus.jl")
include("algorithms.jl")
include("alg_utils.jl")
include("caches/sdirk_caches.jl")
include("caches/implicit_split_step_caches.jl")
include("caches/kencarp_caches.jl")
include("perform_step/sdirk.jl")
include("perform_step/implicit_split_step.jl")
include("perform_step/kencarp.jl")

export ImplicitEM, ImplicitEulerHeun, ImplicitRKMil, STrapezoid, SImplicitMidpoint,
    ISSEM, ISSEulerHeun, SKenCarp

end # module
