module OrdinaryDiffEqPRK

import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache, constvalue, @cache,
    alg_cache, get_fsalfirstlast,
    unwrap_alg, perform_step!, @threaded, isthreaded,
    Sequential, BaseThreads, PolyesterThreads
import SciMLBase: alg_order
import DiffEqBase: initialize!
import MuladdMacro: @muladd
import FastBroadcast: @..

include("algorithms.jl")
include("alg_utils.jl")
include("prk_caches.jl")
include("prk_perform_step.jl")

export KuttaPRK2p5

@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :Sequential, :BaseThreads, :PolyesterThreads))
end

end
