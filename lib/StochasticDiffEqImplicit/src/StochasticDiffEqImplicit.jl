module StochasticDiffEqImplicit

using Reexport: Reexport, @reexport
@reexport using StochasticDiffEqCore
import StochasticDiffEqCore
import DiffEqBase

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, issplit, default_controller, PIController

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step, supports_regular_jumps, isadaptive,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqNewtonAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    StochasticDiffEqCompositeAlgorithm,
    unwrap_alg,
    @cache

import DiffEqBase: @.., SplitSDEFunction, initialize!
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd
import SciMLBase
import SciMLBase: is_diagonal_noise, has_Wfact, _vec, _reshape, _unwrap_val

using OrdinaryDiffEqCore: current_extrapolant!,
    isnewton, set_new_W!, get_W
using OrdinaryDiffEqNonlinearSolve: nlsolvefail, nlsolve!, build_nlsolver,
    markfirststage!, NLNewton
import OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqDifferentiation: calc_J, calc_J!, dolinsolve

using JumpProcesses: JumpProblem
using LinearAlgebra: LinearAlgebra, mul!, norm
using StaticArrays: StaticArrays, SArray
using RecursiveArrayTools: RecursiveArrayTools
using ForwardDiff: ForwardDiff
using FiniteDiff: FiniteDiff
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
