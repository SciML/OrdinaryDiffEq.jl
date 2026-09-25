module OrdinaryDiffEqPDIRK

import OrdinaryDiffEqCore: isfsal,
    OrdinaryDiffEqNewtonAlgorithm, OrdinaryDiffEqConstantCache,
    OrdinaryDiffEqMutableCache, constvalue, alg_cache,
    unwrap_alg, @cache,
    @threaded, perform_step!, isthreaded,
    Sequential, BaseThreads, PolyesterThreads,
    get_fsalfirstlast,
    _fixup_ad
import SciMLBase: alg_order, _unwrap_val
import DiffEqBase: initialize!
import MuladdMacro: @muladd
import FastBroadcast: @..


import SciMLBase

using OrdinaryDiffEqNonlinearSolve: NLNewton, build_nlsolver, nlsolve!, nlsolvefail,
    markfirststage!

import ADTypes: AutoForwardDiff

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
include("pdirk_caches.jl")
include("pdirk_perform_step.jl")

export PDIRK44

@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, :Sequential, :BaseThreads, :PolyesterThreads))
end

end
