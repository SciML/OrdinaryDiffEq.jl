module OrdinaryDiffEqFeagin

import OrdinaryDiffEqCore: perform_step!,
    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqAdaptiveAlgorithm, CompiledFloats,
    alg_cache, @cache,
    constvalue, get_fsalfirstlast,
    trivial_limiter!
import SciMLBase: alg_order
import DiffEqBase: initialize!, calculate_residuals, calculate_residuals!
import FastBroadcast: @..
import MuladdMacro: @muladd
import RecursiveArrayTools: recursivefill!
using DiffEqBase: @tight_loop_macros
import OrdinaryDiffEqCore

using Reexport: Reexport, @reexport
@reexport using SciMLBase: ODEProblem, ODEFunction, solve, init, solve!, step!, remake,
    reinit!, ReturnCode, ContinuousCallback, DiscreteCallback, VectorContinuousCallback,
    CallbackSet, terminate!, add_tstop!, derivative_discontinuity!, set_proposed_dt!,
    successful_retcode, ODEAliasSpecifier
using SciMLBase: SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("feagin_tableaus.jl")
include("feagin_caches.jl")
include("feagin_rk_perform_step.jl")

export Feagin10, Feagin12, Feagin14

end
