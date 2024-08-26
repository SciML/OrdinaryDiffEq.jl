module OrdinaryDiffEqFeagin

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, @unpack, unwrap_alg,
                           calculate_residuals,
                           OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                           constvalue, _unwrap_val, get_fsalfirstlast,
                           generic_solver_docstring, trivial_limiter!,
                           _ode_interpolant!, _ode_addsteps!
using FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def, @tight_loop_macros
using Static: False
import OrdinaryDiffEqCore

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("feagin_tableaus.jl")
include("feagin_caches.jl")
include("feagin_rk_perform_step.jl")

export Feagin10, Feagin12, Feagin14

end
