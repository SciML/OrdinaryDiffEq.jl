module OrdinaryDiffEqRKN

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, alg_extrapolates,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptivePartitionedAlgorithm,
                       OrdinaryDiffEqPartitionedAlgorithm,
                       OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new,
                       trivial_limiter!, _ode_interpolant!, _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def, @tight_loop_macros

include("algorithms.jl")
include("alg_utils.jl")
include("rkn_tableaus.jl")
include("rkn_caches.jl")
include("interp_func.jl")
include("interpolants.jl")
include("rkn_perform_step.jl")

export Nystrom4, FineRKN4, FineRKN5, Nystrom4VelocityIndependent,
       Nystrom5VelocityIndependent,
       IRKN3, IRKN4, DPRKN4, DPRKN5, DPRKN6, DPRKN6FM, DPRKN8, DPRKN12, ERKN4, ERKN5, ERKN7,
       RKN4

end