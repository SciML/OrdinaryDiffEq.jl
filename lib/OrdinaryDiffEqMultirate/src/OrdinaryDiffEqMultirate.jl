module OrdinaryDiffEqMultirate

import OrdinaryDiffEqCore: TmpCache,
    alg_order, isfsal,
    OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    generic_solver_docstring,
    unwrap_alg, initialize!, perform_step!,
    calculate_residuals, calculate_residuals!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    @cache, alg_cache, full_cache, get_fsalfirstlast,
    nlsolve_f, issplit, _fixup_ad, _unwrap_val
import OrdinaryDiffEqCore
import FastBroadcast: @..
import MuladdMacro: @muladd
import RecursiveArrayTools: recursivefill!
import DiffEqBase: prepare_alg
import LinearAlgebra
using SciMLBase: SplitFunction
using OrdinaryDiffEqNonlinearSolve: build_nlsolver, nlsolve!, nlsolvefail,
    markfirststage!, isnewton, set_new_W!, NLNewton
import ADTypes: AutoForwardDiff

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("multirate_tableaus.jl")
include("multirate_caches.jl")
include("multirate_perform_step.jl")

export MREEF, MRAB, MRIGARKERK22a, MRIGARKERK22b, MRIGARKERK33a, MRIGARKERK45a,
    MRIGARKIRK21a, MRIGARKESDIRK34a, MIS

end
