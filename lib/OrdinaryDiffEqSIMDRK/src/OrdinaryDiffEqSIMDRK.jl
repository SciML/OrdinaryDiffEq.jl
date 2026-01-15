module OrdinaryDiffEqSIMDRK

using MuladdMacro, Static
using OrdinaryDiffEqCore: OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqConstantCache,
    trivial_limiter!, calculate_residuals, constvalue
import OrdinaryDiffEqCore: initialize!, perform_step!, alg_cache

using Reexport: @reexport
@reexport using OrdinaryDiffEqCore

include("algorithms.jl")
include("caches.jl")
include("perform_step.jl")

export MER5v2, MER6v2, RK6v4

end
