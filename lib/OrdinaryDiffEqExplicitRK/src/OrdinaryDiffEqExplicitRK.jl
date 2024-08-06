module OrdinaryDiffEqExplicitRK

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order, alg_stability_size, OrdinaryDiffEqAdaptiveAlgorithm,
                       @cache, alg_cache, OrdinaryDiffEqConstantCache, @unpack, unwrap_alg,
                       OrdinaryDiffEqMutableCache, initialize!, perform_step!, isfsal,
                       CompositeAlgorithm, calculate_residuals!, calculate_residuals
using TruncatedStacktraces, RecursiveArrayTools, FastBroadcast, MuladdMacro, DiffEqBase
import LinearAlgebra: norm

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("explicit_rk_caches.jl")
include("explicit_rk_perform_step.jl")

export ExplicitRK

end
