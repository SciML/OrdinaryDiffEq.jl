module OrdinaryDiffEqExponentialRK

import OrdinaryDiffEq: alg_order, alg_adaptive_order, ismultistep, OrdinaryDiffEqExponentialAlgorithm,
                       _unwrap_val, OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                       build_jac_config, UJacobianWrapper, @cache, alg_cache, UDerivativeWrapper,
                       initialize!, perform_step!, @unpack, unwrap_alg, calc_J, calc_J!,
                       OrdinaryDiffEqAdaptiveExponentialAlgorithm, CompositeAlgorithm,
                       ExponentialAlgorithm
using RecursiveArrayTools
using MuladdMacro, FastBroadcast
using LinearAlgebra: axpy!, mul!
using DiffEqBase, SciMLBase
using ExponentialUtilities

include("algorithms.jl")
include("alg_utils.jl")
include("exponential_rk_caches.jl")
include("exponential_rk_perform_step.jl")

export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
       EPIRK4s3B,
       EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43
end
