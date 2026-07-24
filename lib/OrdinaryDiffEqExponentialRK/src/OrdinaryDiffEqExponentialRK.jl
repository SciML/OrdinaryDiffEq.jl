module OrdinaryDiffEqExponentialRK

import OrdinaryDiffEqCore: TmpCache,
    alg_adaptive_order, ismultistep,
    OrdinaryDiffEqExponentialAlgorithm,
    OrdinaryDiffEqMutableCache,
    OrdinaryDiffEqConstantCache,
    @cache, alg_cache,
    perform_step!, unwrap_alg,
    OrdinaryDiffEqAdaptiveExponentialAlgorithm,
    ExponentialAlgorithm, fsal_typeof, isdtchangeable,
    full_cache, get_fsalfirstlast,
    generic_solver_docstring, _fixup_ad, trivial_limiter!
import OrdinaryDiffEqCore
using RecursiveArrayTools: RecursiveArrayTools
import RecursiveArrayTools: recursivecopy!
using MuladdMacro: MuladdMacro, @muladd
using FastBroadcast: FastBroadcast, @..
using LinearAlgebra: axpy!, mul!
import DiffEqBase
import DiffEqBase: calculate_residuals, calculate_residuals!, initialize!
using ExponentialUtilities: ExponentialUtilities, ExpvCache, KrylovSubspace,
    PhivCache, arnoldi, arnoldi!, expv, expv!, phi,
    phiv, phiv!, phiv_timestep, phiv_timestep!
# alg_order / _unwrap_val are owned by SciMLBase (re-exported through OrdinaryDiffEqCore);
# import from the owner to satisfy ExplicitImports' via-owners check.
import SciMLBase: alg_order, _unwrap_val, UJacobianWrapper, UDerivativeWrapper
using OrdinaryDiffEqDifferentiation: build_jac_config, calc_J, calc_J!
import ADTypes: AutoForwardDiff

using Reexport: Reexport, @reexport
@reexport using SciMLBase
using SciMLBase: SciMLBase, SplitFunction

include("algorithms.jl")
include("alg_utils.jl")
include("exponential_rk_caches.jl")
include("exponential_rk_perform_step.jl")

export LawsonEuler, NorsettEuler, ETD1, ETDRK2, ETDRK3, ETDRK4, HochOst4, Exp4, EPIRK4s3A,
    EPIRK4s3B,
    EPIRK5s3, EXPRB53s3, EPIRK5P1, EPIRK5P2, ETD2, Exprb32, Exprb43
end
