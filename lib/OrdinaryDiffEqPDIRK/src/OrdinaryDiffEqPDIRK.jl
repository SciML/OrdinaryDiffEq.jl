module OrdinaryDiffEqPDIRK

import OrdinaryDiffEqCore: isfsal,
    OrdinaryDiffEqNewtonAlgorithm, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqMutableCache, constvalue, alg_cache,
    unwrap_alg, @cache,
    @threaded, perform_step!, isthreaded,
    Sequential, BaseThreads, PolyesterThreads,
    get_fsalfirstlast,
    _fixup_ad
import SciMLBase: alg_order, _unwrap_val
import DiffEqBase: initialize!
import MuladdMacro: @muladd
import FastBroadcast: @..


import SciMLBase

using OrdinaryDiffEqNonlinearSolve: NLNewton, build_nlsolver, nlsolve!, nlsolvefail,
    markfirststage!

import ADTypes: AutoForwardDiff

using Reexport: Reexport, @reexport
@reexport using SciMLBase: ODEProblem, ODEFunction, solve, init, solve!, step!, remake,
    reinit!, ReturnCode, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    CallbackSet, terminate!

include("algorithms.jl")
include("alg_utils.jl")
include("pdirk_caches.jl")
include("pdirk_perform_step.jl")

export PDIRK44

@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :Sequential, :BaseThreads, :PolyesterThreads))
end

end
