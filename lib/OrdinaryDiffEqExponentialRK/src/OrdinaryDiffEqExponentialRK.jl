module OrdinaryDiffEqExponentialRK

import OrdinaryDiffEq: alg_order, alg_adaptive_order, ismultistep, OrdinaryDiffEqExponentialAlgorithm,
                       _unwrap_val, OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                        @cache, alg_cache,
                       initialize!, perform_step!, @unpack, unwrap_alg,
                       OrdinaryDiffEqAdaptiveExponentialAlgorithm, CompositeAlgorithm,
                       ExponentialAlgorithm, fsal_typeof, isdtchangeable, calculate_residuals, calculate_residuals!
using RecursiveArrayTools
using MuladdMacro, FastBroadcast
using LinearAlgebra: axpy!, mul!
using DiffEqBase, SciMLBase
using ExponentialUtilities
import RecursiveArrayTools: recursivecopy!
using OrdinaryDiffEq.OrdinaryDiffEqDifferentiation: build_jac_config, UJacobianWrapper, UDerivativeWrapper, calc_J, calc_J!

include("algorithms.jl")
include("alg_utils.jl")
include("exponential_rk_caches.jl")
include("exponential_rk_perform_step.jl")

export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
       EPIRK4s3B,
       EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43
end
