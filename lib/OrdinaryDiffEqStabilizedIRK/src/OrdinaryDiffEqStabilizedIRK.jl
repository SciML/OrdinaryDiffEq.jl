module OrdinaryDiffEqStabilizedIRK

import OrdinaryDiffEqCore: TmpCache, build_tmp_cache,
    default_controller, IController,
    gamma_default, issplit,
    perform_step!, unwrap_alg, fac_default_gamma,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    alg_cache, @cache,
    get_fsalfirstlast, increment_nf!, set_EEst!,
    generic_solver_docstring, _fixup_ad,
    get_W, isnewton

import SciMLBase: alg_order, _unwrap_val, _reshape, _vec
import DiffEqBase: calculate_residuals, calculate_residuals!, initialize!

import OrdinaryDiffEqDifferentiation: dolinsolve, update_W!
import OrdinaryDiffEqNonlinearSolve: NLNewton, nlsolve!, build_nlsolver,
    markfirststage!, du_alias_or_new

using OrdinaryDiffEqStabilizedRK: ESERK4, ESERK5, SERK2

import FastBroadcast: @..
import RecursiveArrayTools: recursivefill!
import ADTypes: AutoForwardDiff

using Reexport: Reexport, @reexport
using SciMLBase: SciMLBase
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("irkc_utils.jl")
include("irkc_caches.jl")
include("irkc_perform_step.jl")

export IRKC

end
