module OrdinaryDiffEqQPRK

import OrdinaryDiffEq: OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqConstantCache,
                       explicit_rk_docstring, @cache, @unpack, OrdinaryDiffEqMutableCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, @fold, @OnDemandTableauExtract,
                       trivial_limiter!, alg_cache, alg_order, initialize!, perform_step!,
                       constvalue, calculate_residuals!
using Static: False
using MuladdMacro, FastBroadcast
using RecursiveArrayTools: recursive_unitless_bottom_eltype, recursivefill!
include("algorithms.jl")
include("alg_utils.jl")
include("qprk_caches.jl")
include("qprk_tableaus.jl")
include("qprk_perform_step.jl")

export QPRK98

end
