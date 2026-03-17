module StochasticDiffEqWeak

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!, issplit

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, is_split_step, Ihat2,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqNewtonAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    StochasticDiffEqCompositeAlgorithm,
    unwrap_alg,
    calc_twopoint_random!, calc_twopoint_random,
    calc_threepoint_random!, calc_threepoint_random,
    @cache

import DiffEqBase: is_diagonal_noise, @..
import DiffEqBase: calculate_residuals, calculate_residuals!
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd
import SciMLBase

using OrdinaryDiffEqCore: _vec, _reshape, set_new_W!, get_W
using OrdinaryDiffEqNonlinearSolve: NLSolver, nlsolvefail, nlsolve!, build_nlsolver,
    markfirststage!, NLNewton

import OrdinaryDiffEqNonlinearSolve

using LinearAlgebra
using StaticArrays
using RecursiveArrayTools
using ForwardDiff, FiniteDiff
import ADTypes

include("algorithms.jl")
include("alg_utils.jl")

include("caches/srk_weak_caches.jl")
include("caches/implicit_weak_caches.jl")

include("weak_utils.jl")

include("perform_step/srk_weak.jl")
include("perform_step/implicit_weak.jl")

export DRI1, DRI1NM, RI1, RI3, RI5, RI6,
    RDI1WM, RDI2WM, RDI3WM, RDI4WM,
    W2Ito1, RS1, RS2, PL1WM, PL1WMA,
    NON, NON2, COM,
    SIEA, SIEB, SMEA, SMEB,
    IRI1

end # module
