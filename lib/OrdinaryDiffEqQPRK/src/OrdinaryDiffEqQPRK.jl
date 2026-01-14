module OrdinaryDiffEqQPRK

import OrdinaryDiffEqCore: ODEVerbosity, OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqConstantCache,
    explicit_rk_docstring, @cache,
    OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqAdaptiveAlgorithm, @fold, @OnDemandTableauExtract,
    trivial_limiter!, alg_cache, alg_order, initialize!,
    perform_step!, get_fsalfirstlast,
    constvalue, calculate_residuals!, calculate_residuals,
    full_cache
using Static: False
using MuladdMacro, FastBroadcast
using RecursiveArrayTools: recursive_unitless_bottom_eltype, recursivefill!
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("qprk_caches.jl")
include("qprk_tableaus.jl")
include("qprk_perform_step.jl")

export QPRK98

end
