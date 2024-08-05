module OrdinaryDiffEqSSPRK

import OrdinaryDiffEq: alg_order, calculate_residuals!,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       calculate_residuals, ssp_coefficient,
                       OrdinaryDiffEqAlgorithm,
                       OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                       OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
                       OrdinaryDiffEqAdaptiveAlgorithm, uses_uprev,
                       alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                       constvalue, _unwrap_val,
                       explicit_rk_docstring, trivial_limiter!,
                       _ode_interpolant, _ode_interpolant!,
                       _ode_addsteps!
using DiffEqBase, FastBroadcast, Polyester, MuladdMacro, RecursiveArrayTools
using DiffEqBase: @def
using Static: False

import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA

include("algorithms.jl")
include("alg_utils.jl")
include("ssprk_caches.jl")
include("interp_func.jl")
include("ssprk_perform_step.jl")
include("interpolants.jl")
include("addsteps.jl")
include("functions.jl")

export SSPRK53_2N2, SSPRK22, SSPRK53, SSPRK63, SSPRK83, SSPRK43, SSPRK432, SSPRKMSVS32,
       SSPRK54, SSPRK53_2N1, SSPRK104, SSPRK932, SSPRKMSVS43, SSPRK73, SSPRK53_H,
       SSPRK33, SHLDDRK_2N, KYKSSPRK42, SHLDDRK52

end
