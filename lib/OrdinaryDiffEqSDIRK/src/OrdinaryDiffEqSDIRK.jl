module OrdinaryDiffEqSDIRK

import OrdinaryDiffEqCore: alg_order, calculate_residuals!,
                           initialize!, perform_step!, @unpack, unwrap_alg,
                           calculate_residuals, alg_extrapolates,
                           OrdinaryDiffEqAlgorithm,
                           OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                           OrdinaryDiffEqNewtonAdaptiveAlgorithm,
                           OrdinaryDiffEqNewtonAlgorithm,
                           DEFAULT_PRECS,
                           OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats, uses_uprev,
                           alg_cache, _vec, _reshape, @cache, isfsal, full_cache,
                           constvalue, _unwrap_val, _ode_interpolant,
                           trivial_limiter!, _ode_interpolant!,
                           isesdirk, issplit,
                           ssp_coefficient, get_fsalfirstlast, generic_solver_docstring
using TruncatedStacktraces, MuladdMacro, MacroTools, FastBroadcast, RecursiveArrayTools
using SciMLBase: SplitFunction
using LinearAlgebra: mul!, I
import OrdinaryDiffEqCore

using OrdinaryDiffEqDifferentiation: UJacobianWrapper, dolinsolve
using OrdinaryDiffEqNonlinearSolve: du_alias_or_new, markfirststage!, build_nlsolver,
                                    nlsolve!, nlsolvefail, isnewton, get_W, set_new_W!,
                                    NLNewton, COEFFICIENT_MULTISTEP

using Reexport
@reexport using DiffEqBase

include("algorithms.jl")
include("alg_utils.jl")
include("sdirk_caches.jl")
include("kencarp_kvaerno_caches.jl")
include("sdirk_perform_step.jl")
include("kencarp_kvaerno_perform_step.jl")
include("sdirk_tableaus.jl")

export ImplicitEuler, ImplicitMidpoint, Trapezoid, TRBDF2, SDIRK2, SDIRK22,
       Kvaerno3, KenCarp3, Cash4, Hairer4, Hairer42, SSPSDIRK2, Kvaerno4,
       Kvaerno5, KenCarp4, KenCarp47, KenCarp5, KenCarp58, ESDIRK54I8L2SA, SFSDIRK4,
       SFSDIRK5, CFNLIRK3, SFSDIRK6, SFSDIRK7, SFSDIRK8, Kvaerno5, KenCarp4, KenCarp5,
       SFSDIRK4, SFSDIRK5, CFNLIRK3, SFSDIRK6,
       SFSDIRK7, SFSDIRK8, ESDIRK436L2SA2, ESDIRK437L2SA, ESDIRK547L2SA2, ESDIRK659L2SA

end
