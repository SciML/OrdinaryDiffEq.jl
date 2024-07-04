module OrdinaryDiffEqVerner

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new,
                       explicit_rk_docstring, trivial_limiter!,
                       _ode_interpolant!, _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def, @tight_loop_macros
using Static: False

include("algorithms.jl")
include("alg_utils.jl")
include("verner_tableaus.jl")
include("verner_caches.jl")
include("verner_addsteps.jl")
include("interp_func.jl")
include("interpolants.jl")
include("verner_rk_perform_step.jl")

export Vern6, Vern7, Vern8, Vern9

end
