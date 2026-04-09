module OrdinaryDiffEqMultirate

import OrdinaryDiffEqCore: alg_order, isfsal,
    OrdinaryDiffEqAdaptiveAlgorithm,
    generic_solver_docstring,
    unwrap_alg, initialize!, perform_step!,
    calculate_residuals, calculate_residuals!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    @cache, alg_cache, full_cache, get_fsalfirstlast
import OrdinaryDiffEqCore
import FastBroadcast: @..
import MuladdMacro: @muladd
import RecursiveArrayTools: recursivefill!
import DiffEqBase: prepare_alg

using Reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("mreef_caches.jl")
include("mreef_perform_step.jl")

export MREEF

end
