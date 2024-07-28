module OrdinaryDiffEqPRK

import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, alg_order, OrdinaryDiffEqMutableCache,
                       OrdinaryDiffEqConstantCache, constvalue, @unpack, @cache, alg_cache,
                       unwrap_alg, perform_step!, @threaded
import MuladdMacro: @muladd
import FastBroadcast: @..
using Polyester

include("algorithms.jl")
include("alg_utils.jl")
include("prk_caches.jl")
include("prk_perform_step.jl")

export KuttaPRK2p5

end