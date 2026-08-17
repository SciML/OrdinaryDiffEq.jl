module OrdinaryDiffEqQPRK

import OrdinaryDiffEqCore: OrdinaryDiffEqAdaptiveAlgorithm, OrdinaryDiffEqConstantCache,
    @cache,
    OrdinaryDiffEqMutableCache,
    @fold, @OnDemandTableauExtract,
    trivial_limiter!, alg_cache,
    perform_step!, get_fsalfirstlast,
    constvalue
import SciMLBase: alg_order
import DiffEqBase: initialize!, calculate_residuals!, calculate_residuals
using FastBroadcast: FastBroadcast, Serial, @..
using MuladdMacro: MuladdMacro, @muladd
using RecursiveArrayTools: recursive_unitless_bottom_eltype, recursivefill!
import OrdinaryDiffEqCore

using Reexport: Reexport, @reexport
using SciMLBase: ODEFunction, init, solve!, step!, remake, reinit!, ReturnCode,
    ContinuousCallback, DiscreteCallback, VectorContinuousCallback, CallbackSet,
    terminate!, add_tstop!, derivative_discontinuity!, set_proposed_dt!,
    successful_retcode, ODEAliasSpecifier
@reexport using SciMLBase: ODEProblem, solve
using SciMLBase: SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("qprk_caches.jl")
include("qprk_tableaus.jl")
include("qprk_perform_step.jl")

export QPRK98

end
