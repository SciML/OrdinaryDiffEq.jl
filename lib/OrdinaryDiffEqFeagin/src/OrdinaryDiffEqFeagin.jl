module OrdinaryDiffEqFeagin

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, ssp_coefficient,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val, du_alias_or_new,
                       explicit_rk_docstring, trivial_limiter!,
                       _ode_interpolant, _ode_interpolant!,
                       _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def
using Static: False

include("algorithms.jl")
include("alg_utils.jl")
include("feagin_caches.jl")
include("feagin_rk_perform_step.jl")
include("feagin_tableaus.jl")

export Feagin10, Feagin12, Feagin14

end
