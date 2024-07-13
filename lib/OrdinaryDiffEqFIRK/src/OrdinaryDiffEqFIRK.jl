module OrdinaryDiffEqFIRK

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, alg_extrapolates,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptivePartitionedAlgorithm,
                       OrdinaryDiffEqPartitionedAlgorithm,
                       OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new, _ode_interpolant,
                       trivial_limiter!, _ode_interpolant!, _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def, @tight_loop_macros

include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("integrator_interface.jl")
include("firk_tableaus.jl")
include("firk_caches.jl")
include("firk_perform_step.jl")

end
