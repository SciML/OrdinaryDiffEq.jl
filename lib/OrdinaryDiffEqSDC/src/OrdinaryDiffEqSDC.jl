module OrdinaryDiffEqSDC

import OrdinaryDiffEqCore: isfsal,
    OrdinaryDiffEqNewtonAdaptiveAlgorithm,
    alg_adaptive_order,
    generic_solver_docstring,
    unwrap_alg, perform_step!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    @cache, alg_cache, get_fsalfirstlast,
    constvalue, _fixup_ad,
    @threaded, isthreaded
import OrdinaryDiffEqCore
# `alg_order` and `full_cache` are owned by SciMLBase and extended here, so they
# need `import`; `initialize!` is owned by DiffEqBase.
import SciMLBase: alg_order, full_cache
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
import FastBroadcast: @..
import MuladdMacro: @muladd
import LinearAlgebra
import RecursiveArrayTools: recursivefill!
using OrdinaryDiffEqNonlinearSolve: build_nlsolver, nlsolve!, nlsolvefail,
    markfirststage!, NLNewton
import ADTypes: AutoForwardDiff
using EnumX: @enumx

using Reexport: Reexport, @reexport
using SciMLBase: SciMLBase, _unwrap_val
@reexport using SciMLBase

include("sdc_tableaus.jl")
include("min_sr_s_coefficients.jl")
include("algorithms.jl")
include("alg_utils.jl")
include("sdc_caches.jl")
include("sdc_perform_step.jl")

export SDC, SDCNodes, SDCQuadrature, SDCSweeper, SDCStepUpdate

end
