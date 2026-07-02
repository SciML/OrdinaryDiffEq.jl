module StochasticDiffEqCore

using Reexport: Reexport, @reexport
@reexport using DiffEqBase

import ADTypes

import OrdinaryDiffEqCore
# SDE/RODE/Jump algorithm and cache supertypes plus the shared solver-loop
# helpers live in OrdinaryDiffEqCore and are not part of its public API; they
# are dispatched on and extended here (see the EI ignore lists in
# test/qa/qa.jl).
import OrdinaryDiffEqCore: ODEIntegrator,
    StochasticDiffEqAlgorithm, StochasticDiffEqAdaptiveAlgorithm,
    StochasticDiffEqCompositeAlgorithm,
    StochasticDiffEqRODEAlgorithm, StochasticDiffEqRODEAdaptiveAlgorithm,
    StochasticDiffEqRODECompositeAlgorithm,
    StochasticDiffEqNewtonAdaptiveAlgorithm, StochasticDiffEqNewtonAlgorithm,
    StochasticDiffEqJumpAlgorithm, StochasticDiffEqJumpAdaptiveAlgorithm,
    StochasticDiffEqJumpNewtonAdaptiveAlgorithm,
    StochasticDiffEqJumpDiffusionAlgorithm, StochasticDiffEqJumpDiffusionAdaptiveAlgorithm,
    StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm,
    StochasticDiffEqCache, StochasticDiffEqConstantCache, StochasticDiffEqMutableCache,
    beta2_default, beta1_default,
    qmin_default, qmax_default, qsteady_min_default,
    qsteady_max_default, issplit,
    is_composite_algorithm, perform_step!, handle_callback_modifiers!

using RecursiveArrayTools: RecursiveArrayTools, ArrayPartition,
    recursive_bottom_eltype, recursive_unitless_bottom_eltype,
    recursive_unitless_eltype, recursivecopy
using DiffEqNoiseProcess: DiffEqNoiseProcess, CompoundPoissonProcess,
    CompoundPoissonProcess!, NoiseProcess, NoiseTransport, RSWM,
    WienerProcess, WienerProcess!
using ArrayInterface: ArrayInterface
using SimpleNonlinearSolve: SimpleNonlinearSolve, SimpleTrustRegion
using ForwardDiff: ForwardDiff
using StaticArrays: StaticArrays, SArray
using MuladdMacro: MuladdMacro, @muladd
using FiniteDiff: FiniteDiff
using Adapt: Adapt, adapt

import DiffEqBase: ODE_DEFAULT_NORM, ODE_DEFAULT_ISOUTOFDOMAIN,
    ODE_DEFAULT_PROG_MESSAGE, ODE_DEFAULT_UNSTABLE_CHECK

using SciMLOperators: MatrixOperator, WOperator

using Logging: Logging
using SparseArrays: SparseArrays, issparse

# DEVerbosity is owned (and made public) by DiffEqBase but also re-exported
# through OrdinaryDiffEqCore; import it from the owner.
using DiffEqBase: DEVerbosity

using LinearAlgebra: LinearAlgebra, I, mul!
using Random: Random

import ForwardDiff.Dual

import FastPower

import DiffEqBase: step!, initialize!

# Integrator Interface
import DiffEqBase: addat!, full_cache, user_cache, u_cache, du_cache,
    rand_cache, ratenoise_cache,
    resize_non_user_cache!, deleteat_non_user_cache!, addat_non_user_cache!,
    terminate!, get_du, get_dt, get_proposed_dt, set_proposed_dt!,
    savevalues!, add_tstop!, add_saveat!, set_reltol!,
    set_abstol!

using DiffEqBase: @..

import SciMLBase
# AbstractDEAlgorithm/AbstractSDDEIntegrator supertypes are owned and exported
# by SciMLBase; is_diagonal_noise is public SciMLBase API.
import SciMLBase: isadaptive, alg_order, AbstractDEAlgorithm,
    AbstractSDDEIntegrator, is_diagonal_noise
using SciMLBase: CallbackSet, DiscreteProblem, DynamicalSDEFunction,
    NonlinearFunction, NonlinearProblem, ODEAliasSpecifier, SDEProblem,
    get_tmp_cache, isinplace, reinit!
using CommonSolve: init, solve, solve!
using SciMLOperators: SciMLOperators
using DiffEqBase: DiffEqBase

using StochasticDiffEqLevyArea: StochasticDiffEqLevyArea, MaxL2, terms_needed

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

# `@cache` is the cache-definition extension macro that downstream StochasticDiffEq
# solver subpackages import to define their `alg_cache` structs plus the associated
# `full_cache`/`rand_cache`/`ratenoise_cache` accessors (mirrors the public `@cache`
# in OrdinaryDiffEqCore). It is made public (not exported) so ExplicitImports'
# public-API checks recognize it as the supported solver-author surface.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Expr(:public, Symbol("@cache")))
end

end # module
