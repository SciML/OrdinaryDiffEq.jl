module OrdinaryDiffEqFunctionMap

import OrdinaryDiffEqCore: isfsal, beta2_default, beta1_default, OrdinaryDiffEqAlgorithm,
    perform_step!, unwrap_alg,
    OrdinaryDiffEqMutableCache,
    alg_cache, @cache, _ode_addsteps!, _ode_interpolant!,
    _ode_interpolant, get_fsalfirstlast,
    OrdinaryDiffEqConstantCache, dt_required,
    isdiscretecache, isdiscretealg,
    trivial_limiter!
using DiffEqBase: DiffEqBase
import DiffEqBase: initialize!
import RecursiveArrayTools: recursivecopy!
import FastBroadcast: @..
import MuladdMacro: @muladd
import OrdinaryDiffEqCore

import SciMLBase
import SciMLBase: alg_order

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
include("functionmap_caches.jl")
include("interp_func.jl")
include("interpolants.jl")
include("functionmap_perform_step.jl")
include("fixed_timestep_perform_step.jl")

# Default algorithm for DiscreteProblem
function SciMLBase.__solve(
        prob::SciMLBase.DiscreteProblem, ::Nothing, args...;
        kwargs...
    )
    return SciMLBase.__solve(prob, FunctionMap(), args...; kwargs...)
end

function SciMLBase.__init(
        prob::SciMLBase.DiscreteProblem, ::Nothing, args...;
        kwargs...
    )
    return SciMLBase.__init(prob, FunctionMap(), args...; kwargs...)
end

export FunctionMap

end
