module OrdinaryDiffEqVerner

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, alg_stability_size,
                       OrdinaryDiffEqAlgorithm,
                       CompositeAlgorithm, AbstractController, PIDController,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val,
                       explicit_rk_docstring, trivial_limiter!, _ode_interpolant,
                       _ode_interpolant!, _ode_addsteps!, @fold,
                       @OnDemandTableauExtract, AutoAlgSwitch,
                       DerivativeOrderNotPossibleError
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def, @tight_loop_macros
using Static: False
using TruncatedStacktraces
using LinearAlgebra: norm

include("algorithms.jl")
include("alg_utils.jl")
include("verner_tableaus.jl")
include("verner_caches.jl")
include("verner_addsteps.jl")
include("interp_func.jl")
include("interpolants.jl")
include("controllers.jl")
include("verner_rk_perform_step.jl")

export Vern6, Vern7, Vern8, Vern9
export AutoVern6, AutoVern7, AutoVern8, AutoVern9

end
