module StochasticDiffEqCore

using Reexport
@reexport using DiffEqBase

import ADTypes

import OrdinaryDiffEqCore
import OrdinaryDiffEqCore: ODEIntegrator,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCompositeAlgorithm,
    StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
    StochasticDiffEqRODECompositeAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    default_controller, isstandard, ispredictive,
    beta2_default, beta1_default, gamma_default,
    qmin_default, qmax_default, qsteady_min_default,
    qsteady_max_default,
    stepsize_controller!, accept_step_controller,
    step_accept_controller!,
    step_reject_controller!, PIController, DummyController, issplit

# Import shared loop functions from OrdinaryDiffEqCore.
import OrdinaryDiffEqCore: handle_callbacks!, handle_tstop!,
    solution_endpoint_match_cur_integrator!,
    _savevalues!, _postamble!,
    is_composite_cache, is_composite_algorithm, final_progress,
    loopheader!, loopfooter!, _loopfooter!, _step!, perform_step!,
    isaposteriori,
    increment_accept!, increment_reject!,
    calc_dt_propose!, fix_dt_at_bounds!, modify_dt_for_tstops!,
    log_step!, choose_algorithm!, update_uprev!,
    alg_extrapolates, isfsal,
    accept_noise!, reject_noise!, save_noise!, noise_curt, is_noise_saveable,
    reinit_noise!, _determine_initdt, is_constant_cache,
    handle_callback_modifiers!,
    initialize_callbacks!,
    current_extrapolant, current_extrapolant!

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

using OrdinaryDiffEqCore: set_new_W!, get_W, _vec, _reshape

if isdefined(OrdinaryDiffEqCore, :FastConvergence)
    using OrdinaryDiffEqCore:
        FastConvergence, Convergence, SlowConvergence,
        VerySlowConvergence, Divergence

    import OrdinaryDiffEqCore:
        calculate_residuals, calculate_residuals!, nlsolve_f,
        unwrap_cache, islinear
else
    using DiffEqBase:
        FastConvergence, Convergence, SlowConvergence, VerySlowConvergence,
        Divergence

    import DiffEqBase:
        calculate_residuals, calculate_residuals!, nlsolve_f, unwrap_cache,
        islinear
end

import SciMLBase
import SciMLBase: isadaptive, alg_order

using StochasticDiffEqLevyArea

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
include("caches/cache_types.jl")
include("integrators/type.jl")
include("alg_utils.jl")
include("integrators/stepsize_controllers.jl")
include("integrators/integrator_utils.jl")
include("integrators/integrator_interface.jl")
include("initialize_dae.jl")
include("solve.jl")
include("perform_step/composite.jl")
include("iterated_integrals.jl")
include("composite_algs.jl")
include("weak_utils.jl")

export StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCompositeAlgorithm,
    StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
    StochasticDiffEqRODECompositeAlgorithm

export StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqNewtonAlgorithm,
    StochasticDiffEqJumpAlgorithm, StochasticDiffEqJumpAdaptiveAlgorithm,
    StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
    StochasticDiffEqJumpDiffusionAlgorithm, StochasticDiffEqJumpDiffusionAdaptiveAlgorithm,
    StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm

export IteratedIntegralApprox, IICommutative, IILevyArea

export StochasticCompositeAlgorithm, AutoSwitch, AutoAlgSwitch

export SDEIntegrator, SDEOptions, SDEAlgTypes

# Re-export cache types for solver subpackages
export StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    StochasticCompositeCache

# Export core utility functions used by solver subpackages
export alg_cache, alg_order, get_current_alg_order, isadaptive,
    alg_compatible, alg_needs_extra_process, supports_regular_jumps,
    alg_mass_matrix_compatible, alg_can_repeat_jac,
    is_split_step, alg_control_rate, delta_default,
    alg_stability_size, unwrap_alg, issplit,
    DiffCache, get_du, get_chunksize, determine_chunksize,
    NLSOLVEJL_SETUP, IIFNLSolveFunc, DiffEqNLSolveTag

# Export iterated integral utilities for solver subpackages
export get_iterated_I, get_iterated_I!, get_Jalg,
    AbstractJ, AbstractJDiagonal, AbstractJCommute,
    JDiagonal_oop, JDiagonal_iip, JCommute_oop, JCommute_iip,
    IteratedIntegralAlgorithm_iip

# Export weak utility functions
export calc_twopoint_random!, calc_twopoint_random,
    calc_threepoint_random!, calc_threepoint_random, Ihat2

# Export noise interface
export resize_noise!, fill_new_noise_caches!,
    deleteat_noise!, addat_noise!

# Export Z prototype dispatch for noise creation
export _z_prototype, concrete_prob, _resolve_rng, _sde_init

# Export TauLeapingDrift for Jump subpackage
export TauLeapingDrift

# General functions
export solve, init, solve!, step!

# Export misc tools (check functions)
export checkSRIOrder, checkSRAOrder

end # module
