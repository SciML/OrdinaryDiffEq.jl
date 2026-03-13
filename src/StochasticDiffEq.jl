"""
$(DocStringExtensions.README)
"""
module StochasticDiffEq

using DocStringExtensions

using Reexport
@reexport using DiffEqBase

import ADTypes

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: ODEIntegrator,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCompositeAlgorithm,
    StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
    StochasticDiffEqRODECompositeAlgorithm,
    default_controller, isstandard, ispredictive,
    beta2_default, beta1_default, gamma_default,
    qmin_default, qmax_default, qsteady_min_default,
    qsteady_max_default,
    stepsize_controller!, accept_step_controller,
    step_accept_controller!,
    step_reject_controller!, PIController, DummyController, issplit

# Import shared loop functions from OrdinaryDiffEqCore.
# SDE reuses ODE's unified loop functions directly — no hook overrides needed.
# SDE-specific behavior is handled via field conditionals (isnothing(integrator.W))
# and algorithm traits (isfsal, alg_extrapolates, etc.).
import OrdinaryDiffEqCore: handle_callbacks!, handle_tstop!,
    solution_endpoint_match_cur_integrator!,
    _savevalues!, _postamble!,
    is_composite_cache, is_composite_algorithm, final_progress,
    loopheader!, loopfooter!, _loopfooter!, _step!, perform_step!,
    isaposteriori, fixed_t_for_floatingpoint_error!,
    increment_accept!, increment_reject!,
    calc_dt_propose!, fix_dt_at_bounds!, modify_dt_for_tstops!,
    log_step!, choose_algorithm!, update_uprev!,
    alg_extrapolates, isfsal,
    accept_noise!, reject_noise!, save_noise!, noise_curt, is_noise_saveable,
    handle_callback_modifiers!,
    initialize_callbacks!

using RecursiveArrayTools, DataStructures
using DiffEqNoiseProcess, Random, ArrayInterface
using SimpleNonlinearSolve, ForwardDiff, StaticArrays, MuladdMacro, FiniteDiff, Base.Threads
using Adapt

import DiffEqBase: ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN,
    ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

using SciMLOperators: MatrixOperator, WOperator

using DiffEqBase: TimeGradientWrapper, UJacobianWrapper, TimeDerivativeWrapper,
    UDerivativeWrapper

import RecursiveArrayTools: chain

using Logging, SparseArrays

using SciMLLogging: AbstractVerbositySpecifier, AbstractVerbosityPreset,
    AbstractMessageLevel, None, Minimal, Standard, Detailed, All,
    Silent, DebugLevel, InfoLevel, WarnLevel, ErrorLevel, @SciMLMessage

using OrdinaryDiffEqCore: DEVerbosity

using LinearAlgebra, Random

import ForwardDiff.Dual

import FastPower

import DiffEqBase: step!, initialize!, DEAlgorithm,
    AbstractSDEAlgorithm, AbstractRODEAlgorithm, DEIntegrator,
    AbstractDiffEqInterpolation,
    DECache, AbstractSDEIntegrator, AbstractRODEIntegrator,
    AbstractContinuousCallback,
    Tableau, AbstractSDDEIntegrator

# Integrator Interface
import DiffEqBase: resize!, deleteat!, addat!, full_cache, user_cache, u_cache, du_cache,
    rand_cache, ratenoise_cache,
    resize_non_user_cache!, deleteat_non_user_cache!, addat_non_user_cache!,
    terminate!, get_du, get_dt, get_proposed_dt, set_proposed_dt!,
    u_modified!, savevalues!, add_tstop!, add_saveat!, set_reltol!,
    set_abstol!, postamble!, last_step_failed, has_Wfact, has_jac,
    get_tstops, get_tstops_array, get_tstops_max

using DiffEqBase: check_error!, is_diagonal_noise, @..

using OrdinaryDiffEqNonlinearSolve: nlsolvefail, isnewton
using OrdinaryDiffEqCore: set_new_W!, get_W, _vec, _reshape

using OrdinaryDiffEqNonlinearSolve: NLSolver

import OrdinaryDiffEqCore
import OrdinaryDiffEqDifferentiation
import OrdinaryDiffEqNonlinearSolve

if isdefined(OrdinaryDiffEqCore, :FastConvergence)
    using OrdinaryDiffEqCore:
        FastConvergence, Convergence, SlowConvergence,
        VerySlowConvergence, Divergence

    import OrdinaryDiffEqCore:
        calculate_residuals, calculate_residuals!, nlsolve_f,
        unwrap_cache, islinear

    using OrdinaryDiffEqNonlinearSolve: NLFunctional, NLAnderson, NLNewton
else
    using DiffEqBase:
        FastConvergence, Convergence, SlowConvergence, VerySlowConvergence,
        Divergence

    import DiffEqBase:
        calculate_residuals, calculate_residuals!, nlsolve_f, unwrap_cache,
        islinear
end

import SciMLBase

using LevyArea

const CompiledFloats = Union{Float32, Float64}

import JumpProcesses
import JumpProcesses: JumpProblem

import Base.Threads
@static if VERSION < v"1.3"
    seed_multiplier() = Threads.threadid()
else
    seed_multiplier() = 1
end

include("misc_utils.jl")
include("algorithms.jl")
include("interp_func.jl")
include("caches/cache_types.jl")
include("caches/basic_method_caches.jl")
include("caches/explicit_3s_mil_methods.jl")
include("caches/lamba_caches.jl")
include("caches/iif_caches.jl")
include("caches/sdirk_caches.jl")
include("caches/implicit_split_step_caches.jl")
include("caches/sra_caches.jl")
include("caches/rossler_caches.jl")
include("caches/srk_weak_caches.jl")
include("caches/kencarp_caches.jl")
include("caches/implicit_weak_caches.jl")
include("caches/predcorr_caches.jl")
include("caches/SROCK_caches.jl")
include("caches/tau_caches.jl")
include("caches/dynamical_caches.jl")
include("integrators/type.jl")
include("dense.jl")
include("alg_utils.jl")
include("integrators/stepsize_controllers.jl")
include("integrators/integrator_utils.jl")
include("cache_utils.jl")
include("integrators/integrator_interface.jl")
include("iterator_interface.jl")
include("initialize_dae.jl")
include("solve.jl")
include("initdt.jl")
include("perform_step/low_order.jl")
include("perform_step/explicit_3s_mil_methods.jl")
include("perform_step/lamba.jl")
include("perform_step/iif.jl")
include("perform_step/sri.jl")
include("perform_step/sra.jl")
include("perform_step/srk_weak.jl")
include("perform_step/sdirk.jl")
include("perform_step/implicit_split_step.jl")
include("perform_step/kencarp.jl")
include("perform_step/implicit_weak.jl")
include("perform_step/predcorr.jl")
include("perform_step/split.jl")
include("perform_step/composite.jl")
include("perform_step/SROCK_perform_step.jl")
include("perform_step/tau_leaping.jl")
include("perform_step/dynamical.jl")
include("tableaus.jl")
include("SROCK_tableaus.jl")
include("iterated_integrals.jl")
include("SROCK_utils.jl")
include("composite_algs.jl")
include("weak_utils.jl")
include("default_sde_alg.jl")

export StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticCompositeAlgorithm

export EM, LambaEM, PCEuler, RKMil, SRA, SRI, SRIW1,
    SRA1, SRA2, SRA3,
    SOSRA, SOSRA2, RKMilCommute, RKMilGeneral,
    SRIW2, SOSRI, SOSRI2, SKenCarp,
    SROCK1, SROCK2, SROCKEM, SKSROCK, TangXiaoSROCK2, KomBurSROCK2, SROCKC2,
    WangLi3SMil_A, WangLi3SMil_B, WangLi3SMil_C, WangLi3SMil_D, WangLi3SMil_E,
    WangLi3SMil_F,
    AutoSOSRI2, AutoSOSRA2,
    DRI1, DRI1NM, RI1, RI3, RI5, RI6, RDI1WM, RDI2WM, RDI3WM, RDI4WM, W2Ito1,
    RS1, RS2,
    PL1WM, PL1WMA,
    NON, COM, NON2

export SIEA, SMEA, SIEB, SMEB

export EulerHeun, LambaEulerHeun

export SimplifiedEM

export SplitEM, IIF1M, IIF2M, IIF1Mil

export ImplicitEM, ImplicitEulerHeun, ISSEM, ISSEulerHeun,
    ImplicitRKMil, STrapezoid, SImplicitMidpoint, IRI1

export TauLeaping, CaoTauLeaping, ImplicitTauLeaping, ThetaTrapezoidalTauLeaping

export BAOAB

export StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
    StochasticDiffEqRODECompositeAlgorithm

export RandomEM, RandomTamedEM, RandomHeun

export IteratedIntegralApprox, IICommutative, IILevyArea

#General Functions
export solve, init, solve!, step!

#Misc Tools
export checkSRIOrder, checkSRAOrder, checkRIOrder, checkRSOrder,
    checkNONOrder,
    constructSRIW1, constructSRA1,
    constructDRI1, constructRI1, constructRI3, constructRI5, constructRI6,
    constructRDI1WM, constructRDI2WM, constructRDI3WM, constructRDI4WM,
    constructRS1, constructRS2,
    constructNON, constructNON2

end # module
