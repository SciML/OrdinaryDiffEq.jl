module OrdinaryDiffEqPRK

import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, alg_order, OrdinaryDiffEqMutableCache,
                           OrdinaryDiffEqConstantCache, constvalue, @unpack, @cache,
                           alg_cache, get_fsalfirstlast,
                           unwrap_alg, perform_step!, @threaded, initialize!, isthreaded,
                           full_cache
import MuladdMacro: @muladd
import FastBroadcast: @..
using Polyester

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("prk_caches.jl")
include("prk_perform_step.jl")

export KuttaPRK2p5

end
