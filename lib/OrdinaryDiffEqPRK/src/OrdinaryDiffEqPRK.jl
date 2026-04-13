module OrdinaryDiffEqPRK

import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, alg_order, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache, constvalue, @cache,
    alg_cache, get_fsalfirstlast,
    unwrap_alg, perform_step!, @threaded, initialize!, isthreaded,
    full_cache, generic_solver_docstring
import MuladdMacro: @muladd
import FastBroadcast: @..
using Polyester

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("prk_caches.jl")
include("prk_perform_step.jl")

export KuttaPRK2p5

end
