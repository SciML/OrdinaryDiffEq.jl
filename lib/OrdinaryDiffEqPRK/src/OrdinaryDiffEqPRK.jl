module OrdinaryDiffEqPRK

import OrdinaryDiffEqCore: TmpCache, build_tmp_cache,
    OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache, constvalue, @cache,
    alg_cache, get_fsalfirstlast,
    unwrap_alg, perform_step!, @threaded, isthreaded,
    generic_solver_docstring
import SciMLBase: alg_order, full_cache
import DiffEqBase: initialize!
import MuladdMacro: @muladd
import FastBroadcast: @..

using Reexport: Reexport, @reexport
using SciMLBase: SciMLBase
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("prk_caches.jl")
include("prk_perform_step.jl")

export KuttaPRK2p5

end
