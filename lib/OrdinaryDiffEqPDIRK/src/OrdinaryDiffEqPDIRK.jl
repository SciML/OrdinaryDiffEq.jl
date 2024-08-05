module OrdinaryDiffEqPDIRK

import OrdinaryDiffEq: isfsal, alg_order, _unwrap_val,
                       OrdinaryDiffEqNewtonAlgorithm, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqMutableCache, constvalue, alg_cache,
                       uses_uprev, @unpack, unwrap_alg, @cache, DEFAULT_PRECS,
                       @threaded, initialize!, perform_step!, isthreaded
import StaticArrays: SVector
import MuladdMacro: @muladd
import FastBroadcast: @..
using Polyester

using OrdinaryDiffEq.OrdinaryDiffEqDifferentiation: dolinsolve
using OrdinaryDiffEq.OrdinaryDiffEqNonlinearSolve: NLNewton, build_nlsolver, nlsolve!, nlsolvefail, markfirststage!

include("algorithms.jl")
include("alg_utils.jl")
include("pdirk_caches.jl")
include("pdirk_perform_step.jl")

export PDIRK44

end
