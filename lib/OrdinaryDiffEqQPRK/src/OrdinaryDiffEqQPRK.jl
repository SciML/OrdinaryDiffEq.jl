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
# Re-exported SciMLBase API. This list is identical in every OrdinaryDiffEq solver
# sublibrary and is generated from the rule documented in docs/src/api/reexports.md;
# `test/qa/qa_tests.jl` checks the sublibraries against each other and against that rule.
@reexport using SciMLBase: AnalyticalProblem, BVProblem, ConvexOptimizationProblem,
    DAEProblem, DDEProblem, DiscreteProblem, DynamicalDDEProblem, DynamicalODEProblem,
    DynamicalSDEProblem, EigenvalueProblem, EnsembleProblem, HomotopyProblem,
    ImmutableNonlinearProblem, ImmutableODEProblem, ImplicitDiscreteProblem,
    IncrementingODEProblem, IntegralProblem, IntervalNonlinearProblem, LinearProblem,
    NoiseProblem, NonlinearLeastSquaresProblem, NonlinearProblem, ODEProblem,
    OptimizationProblem, PDEProblem, RODEProblem, SCCNonlinearProblem, SDDEProblem,
    SDEProblem, SampledIntegralProblem, SecondOrderBVProblem, SecondOrderDDEProblem,
    SecondOrderODEProblem, SplitODEProblem, SplitSDEProblem, SteadyStateProblem,
    TwoPointBVProblem, TwoPointSecondOrderBVProblem, WeightedEnsembleProblem, BVPFunction,
    BatchIntegralFunction, DAEFunction, DDEFunction, DiscreteFunction, DynamicalBVPFunction,
    DynamicalDDEFunction, DynamicalODEFunction, DynamicalSDEFunction,
    HomotopyNonlinearFunction, ImplicitDiscreteFunction, IncrementingODEFunction,
    IntegralFunction, IntervalNonlinearFunction, MultiObjectiveOptimizationFunction,
    NonlinearFunction, ODEFunction, ODEInputFunction, OptimizationFunction, RODEFunction,
    SDDEFunction, SDEFunction, SplitFunction, SplitSDEFunction, TwoPointBVPFunction,
    TwoPointDynamicalBVPFunction, solve, solve!, init, step!, remake, ReturnCode,
    successful_retcode, DEStats, NLStats, NullParameters, AutoSpecialize, FullSpecialize,
    NoSpecialize, FunctionWrapperSpecialize, CheckInit, NoInit, OverrideInit, CallbackSet,
    ContinuousCallback, DiscreteCallback, VectorContinuousCallback, LeftRootFind,
    RightRootFind, NoRootFind, add_saveat!, add_tstop!, auto_dt_reset!,
    change_t_via_interpolation!, check_error, check_error!, first_tstop, get_dt, get_du,
    get_du!, get_proposed_dt, get_tmp_cache, pop_tstop!, reinit!, savevalues!, set_abstol!,
    set_proposed_dt!, set_reltol!, set_t!, set_u!, set_ut!, terminate!, u_modified!,
    EnsembleAnalysis, EnsembleDistributed, EnsembleSerial, EnsembleSplitThreads,
    EnsembleSummary, EnsembleThreads
using SciMLBase: SciMLBase
using SciMLBase: SciMLBase

include("algorithms.jl")
include("alg_utils.jl")
include("qprk_caches.jl")
include("qprk_tableaus.jl")
include("qprk_perform_step.jl")

export QPRK98

end
