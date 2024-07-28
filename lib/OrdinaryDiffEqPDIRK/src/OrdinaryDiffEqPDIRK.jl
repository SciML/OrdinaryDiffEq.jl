module OrdinaryDiffEqPDIRK

import OrdinaryDiffEq: isfsal, alg_order, _unwrap_val, DEFAULT_PRECS, NLNewton,
                       OrdinaryDiffEqNewtonAlgorithm, dolinsolve, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqMutableCache, constvalue, alg_cache, build_nlsolver,
                       uses_uprev, nlsolve!, nlsolvefail, @unpack, unwrap_alg, @cache,
                       markfirststage!, @threaded
import StaticArrays: SVector
import MuladdMacro: @muladd
import FastBroadcast: @..
using Polyester

include("algorithms.jl")
include("alg_utils.jl")
include("pdirk_caches.jl")
include("pdirk_perform_step.jl")

export PDIRK44

end