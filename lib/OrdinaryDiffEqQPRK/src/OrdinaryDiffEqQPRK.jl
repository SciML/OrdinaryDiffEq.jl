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
# Kept in sync with docs/src/api/reexports.md by test/qa/reexport_tests.jl.
@reexport using SciMLBase: DAEProblem, DiscreteProblem, DynamicalODEProblem,
    EnsembleProblem, ODEProblem, SecondOrderODEProblem, SplitODEProblem, DAEFunction,
    DiscreteFunction, DynamicalODEFunction, ODEFunction, SplitFunction, solve, solve!, init,
    step!, remake, ReturnCode, successful_retcode, DEStats, NLStats, NullParameters,
    AutoSpecialize, FullSpecialize, NoSpecialize, FunctionWrapperSpecialize, CheckInit,
    NoInit, OverrideInit, CallbackSet, ContinuousCallback, DiscreteCallback,
    VectorContinuousCallback, LeftRootFind, RightRootFind, NoRootFind, add_saveat!,
    add_tstop!, auto_dt_reset!, change_t_via_interpolation!, check_error, check_error!,
    first_tstop, get_dt, get_du, get_du!, get_proposed_dt, get_tmp_cache, pop_tstop!,
    reinit!, savevalues!, set_abstol!, set_proposed_dt!, set_reltol!, set_t!, set_u!,
    set_ut!, terminate!, u_modified!, EnsembleAnalysis, EnsembleDistributed, EnsembleSerial,
    EnsembleSplitThreads, EnsembleSummary, EnsembleThreads
using SciMLBase: SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("qprk_caches.jl")
include("qprk_tableaus.jl")
include("qprk_perform_step.jl")

export QPRK98

end
