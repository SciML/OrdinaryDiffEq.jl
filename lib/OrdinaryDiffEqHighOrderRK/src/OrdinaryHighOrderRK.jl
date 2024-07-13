module OrdinaryHighOrderRK

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqPartitionedAlgorithm,
                       CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new,
                       explicit_rk_docstring, trivial_limiter!,
                       _ode_interpolant!, _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools

include("algorithms.jl")
include("alg_utils.jl")
include("high_order_rk_caches.jl")
include("interp_func.jl")
include("interpolants.jl")
include("high_order_rk_addsteps.jl")
include("high_order_rk_tableaus.jl")
include("high_order_rk_perform_step.jl")

export TanYam7, DP8, TsitPap8, PFRK87

end 
