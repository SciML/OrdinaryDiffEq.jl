module OrdinaryDiffEqExponentialRK

import OrdinaryDiffEqCore: ODEVerbosity, alg_order, alg_adaptive_order, ismultistep,
    OrdinaryDiffEqExponentialAlgorithm,
    _unwrap_val, OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache,
    @cache, alg_cache,
    initialize!, perform_step!, unwrap_alg,
    OrdinaryDiffEqAdaptiveExponentialAlgorithm, CompositeAlgorithm,
    ExponentialAlgorithm, fsal_typeof, isdtchangeable,
    calculate_residuals, calculate_residuals!,
    full_cache, get_fsalfirstlast,
    generic_solver_docstring, _bool_to_ADType, _process_AD_choice
import OrdinaryDiffEqCore
using RecursiveArrayTools
using MuladdMacro, FastBroadcast
using LinearAlgebra: axpy!, mul!
import DiffEqBase
import DiffEqBase: prepare_alg
using ExponentialUtilities
import RecursiveArrayTools: recursivecopy!
using OrdinaryDiffEqDifferentiation: build_jac_config, UJacobianWrapper, UDerivativeWrapper,
    calc_J, calc_J!
import ADTypes: AutoForwardDiff, AbstractADType

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("exponential_rk_caches.jl")
include("exponential_rk_perform_step.jl")

export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
    EPIRK4s3B,
    EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43
end
