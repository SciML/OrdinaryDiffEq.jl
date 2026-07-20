module OrdinaryDiffEqPDIRK

import OrdinaryDiffEqCore: TmpCache, 
    isfsal,
    OrdinaryDiffEqNewtonAlgorithm, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqMutableCache, constvalue, alg_cache,
    unwrap_alg, @cache,
    @threaded, perform_step!, isthreaded,
    full_cache, get_fsalfirstlast, differentiation_rk_docstring,
    _fixup_ad
import SciMLBase: alg_order, _unwrap_val
import DiffEqBase: initialize!
import MuladdMacro: @muladd
import FastBroadcast: @..


using Reexport: Reexport, @reexport
using SciMLBase: SciMLBase
@reexport using SciMLBase

using OrdinaryDiffEqNonlinearSolve: NLNewton, build_nlsolver, nlsolve!, nlsolvefail,
    markfirststage!

import ADTypes: AutoForwardDiff

include("algorithms.jl")
include("alg_utils.jl")
include("pdirk_caches.jl")
include("pdirk_perform_step.jl")

export PDIRK44

end
