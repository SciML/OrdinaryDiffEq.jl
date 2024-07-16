module OrdinaryDiffEqFIRK

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

include("algorithms.jl")
include("alg_utils.jl")
include("controllers.jl")
include("firk_caches.jl")
include("firk_tableaus.jl")
include("firk_perform_step.jl")
include("integrator_interface.jl")

export RadauIIA3, RadauIIA5

end
