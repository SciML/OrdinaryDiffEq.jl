module OrdinaryDiffEqFeagin

import OrdinaryDiffEqCore: perform_step!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats,
    alg_cache, @cache, full_cache,
    constvalue, get_fsalfirstlast,
    generic_solver_docstring, trivial_limiter!
import SciMLBase: alg_order
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
import FastBroadcast: @..
import MuladdMacro: @muladd
import RecursiveArrayTools: recursivefill!
using DiffEqBase: @tight_loop_macros
import OrdinaryDiffEqCore

using Reexport: @reexport
@reexport using SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("feagin_tableaus.jl")
include("feagin_caches.jl")
include("feagin_rk_perform_step.jl")

export Feagin10, Feagin12, Feagin14

end
