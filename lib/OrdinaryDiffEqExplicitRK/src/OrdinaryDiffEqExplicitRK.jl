module OrdinaryDiffEqExplicitRK

import OrdinaryDiffEqCore: alg_order, alg_adaptive_order, alg_stability_size,
                           OrdinaryDiffEqAdaptiveAlgorithm,
                           @cache, alg_cache, OrdinaryDiffEqConstantCache, @unpack,
                           unwrap_alg,
                           OrdinaryDiffEqMutableCache, initialize!, perform_step!, isfsal,
                           CompositeAlgorithm, calculate_residuals!, calculate_residuals,
                           full_cache, get_fsalfirstlast,
                           _ode_interpolant, _ode_interpolant!,
                           DerivativeOrderNotPossibleError
using TruncatedStacktraces: @truncate_stacktrace
using RecursiveArrayTools, FastBroadcast, MuladdMacro, DiffEqBase
import LinearAlgebra: norm
import OrdinaryDiffEqCore

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("explicit_rk_caches.jl")
include("explicit_rk_perform_step.jl")
include("interp_func.jl")
include("interpolants.jl")

export ExplicitRK

end
