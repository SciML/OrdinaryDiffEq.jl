module StochasticDiffEqLeaping

using Reexport
@reexport using StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, initialize!, issplit,
    stepsize_controller!, accept_step_controller, step_accept_controller!,
    step_reject_controller!, DummyController

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, alg_control_rate,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqJumpAlgorithm, StochasticDiffEqJumpAdaptiveAlgorithm,
    StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    TauLeapingDrift,
    @cache

import DiffEqBase: is_diagonal_noise, @..
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd
import SciMLBase

using OrdinaryDiffEqCore: _vec, _reshape
using OrdinaryDiffEqNonlinearSolve: NLSolver, NLFunctional, NLAnderson, NLNewton,
    nlsolvefail, isnewton

using LinearAlgebra
using StaticArrays
using RecursiveArrayTools

import JumpProcesses

import OrdinaryDiffEqNonlinearSolve

include("algorithms.jl")
include("alg_utils.jl")
include("stepsize_controllers.jl")
include("integrator_utils.jl")
include("caches/tau_caches.jl")
include("perform_step/tau_leaping.jl")

export TauLeaping, CaoTauLeaping, ImplicitTauLeaping, ThetaTrapezoidalTauLeaping

end # module
