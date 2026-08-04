module OrdinaryDiffEqMultirate

import OrdinaryDiffEqCore: alg_order, isfsal,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    OrdinaryDiffEqRosenbrockAdaptiveAlgorithm,
    generic_solver_docstring, differentiation_rk_docstring,
    unwrap_alg, initialize!, perform_step!,
    calculate_residuals, calculate_residuals!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    @cache, alg_cache, get_fsalfirstlast,
    nlsolve_f, issplit, _fixup_ad, _unwrap_val,
    _vec, _reshape, LinearAliasSpecifier
import OrdinaryDiffEqCore
import FastBroadcast: @..
import MuladdMacro: @muladd
import RecursiveArrayTools: recursivefill!
import DiffEqBase
import DiffEqBase: prepare_alg
import LinearAlgebra
import SciMLBase: full_cache
import LinearAlgebra: UniformScaling
import LinearSolve
import SparseArrays: SparseMatrixCSC, getcolptr, issparse, nonzeros, rowvals
using SciMLBase: SplitFunction, LinearProblem, init
import SciMLBase
import SciMLOperators: AbstractSciMLOperator, WOperator, update_coefficients!
using OrdinaryDiffEqNonlinearSolve: build_nlsolver, nlsolve!, nlsolvefail,
    markfirststage!, isnewton, set_new_W!, NLNewton
using OrdinaryDiffEqDifferentiation: build_uf, build_jac_config, build_J_W,
    calc_J, calc_J!, jacobian2W, jacobian2W!, dolinsolve, wrapprecs, issuccess_W
import ADTypes: AutoForwardDiff

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("multirate_tableaus.jl")
include("multirate_caches.jl")
include("multirate_perform_step.jl")

export MREEF, MREIL, MRAB, MRIGARKERK22a, MRIGARKERK22b, MRIGARKERK33a, MRIGARKERK45a,
    MRIGARKIRK21a, MRIGARKESDIRK34a, MRIGARKESDIRK46a, MIS

end
