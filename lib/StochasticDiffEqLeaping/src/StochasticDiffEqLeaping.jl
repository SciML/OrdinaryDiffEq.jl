module StochasticDiffEqLeaping

using Reexport: Reexport, @reexport
@reexport using StochasticDiffEqCore
using StochasticDiffEqCore: StochasticDiffEqCore

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: perform_step!, issplit,
    stepsize_controller!, step_accept_controller!,
    step_reject_controller!, DummyController

import StochasticDiffEqCore: alg_cache, alg_order, alg_compatible,
    alg_needs_extra_process, alg_control_rate,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqJumpAlgorithm, StochasticDiffEqJumpAdaptiveAlgorithm,
    StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    TauLeapingDrift,
    @cache

import DiffEqBase
# `@..` is the SciML fused-broadcast macro; its owner is FastBroadcast (not a
# direct dependency here), so it is re-exported through DiffEqBase.
import DiffEqBase: @..
import DiffEqBase: full_cache, rand_cache, ratenoise_cache

import MuladdMacro: @muladd
import SciMLBase

using OrdinaryDiffEqNonlinearSolve: NLFunctional

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
