module OrdinaryDiffEqCore

if isdefined(Base, :Experimental) &&
        isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

import DocStringExtensions
import Reexport: @reexport
using Reexport: @reexport
@reexport using SciMLBase
import DiffEqBase

import Logging: @logmsg, LogLevel

using MuladdMacro: @muladd

using LinearAlgebra: opnorm, I, UniformScaling, diag, rank, isdiag

import PrecompileTools

import FillArrays: Trues, Falses

import FastPower: fastpower

# Interfaces
import SciMLBase: solve!, step!, isadaptive
import DiffEqBase: initialize!

# DAE Initialization algorithms
import DiffEqBase: DefaultInit, ShampineCollocationInit, BrownFullBasicInit

# Internal utils
import DiffEqBase: ODE_DEFAULT_NORM,
    ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE,
    ODE_DEFAULT_UNSTABLE_CHECK

import SciMLOperators: AbstractSciMLOperator, AbstractSciMLScalarOperator,
    MatrixOperator, FunctionOperator,
    update_coefficients, update_coefficients!, DEFAULT_UPDATE_FUNC,
    isconstant

using DiffEqBase: DEIntegrator

import Random

import RecursiveArrayTools: chain, recursivecopy!, recursivecopy, recursive_bottom_eltype, recursive_unitless_bottom_eltype, recursive_unitless_eltype, copyat_or_push!, DiffEqArray, recursivefill!

import RecursiveArrayTools
using DataStructures: BinaryHeap, FasterForward
import DataStructures
using ArrayInterface: ArrayInterface, issingular

import TruncatedStacktraces: @truncate_stacktrace, VERBOSE_MSG

import StaticArraysCore: SArray, MVector, SVector, StaticArray, MMatrix,
    StaticMatrix

# Integrator Interface
import SciMLBase: resize!, deleteat!, addat!, full_cache, user_cache, u_cache, du_cache,
    resize_non_user_cache!, deleteat_non_user_cache!, addat_non_user_cache!,
    terminate!, get_du, get_dt, get_proposed_dt, set_proposed_dt!,
    u_modified!, savevalues!,
    add_tstop!, has_tstop, first_tstop, pop_tstop!,
    add_saveat!, set_reltol!,
    set_abstol!, postamble!, last_step_failed,
    isautodifferentiable
import DiffEqBase: get_tstops, get_tstops_array, get_tstops_max

using DiffEqBase: check_error!, @def, _vec, _reshape

using FastBroadcast: @.., True, False

using SciMLBase: NoInit, CheckInit, OverrideInit, AbstractDEProblem, _unwrap_val,
    ODEAliasSpecifier

import SciMLBase: AbstractNonlinearProblem, alg_order, LinearAliasSpecifier

import SciMLBase: unwrap_cache,
    islinear
import DiffEqBase: calculate_residuals,
    calculate_residuals!, @tight_loop_macros,
    timedepentdtmin

import Polyester
# MacroTools and Adapt imported but not directly used in OrdinaryDiffEqCore
# using MacroTools, Adapt
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AbstractADType, AutoSparse, dense_ad
import Accessors: @reset

# SciMLStructures symbols imported but not directly used in OrdinaryDiffEqCore
# using SciMLStructures: canonicalize, Tunable, isscimlstructure

using SciMLLogging: SciMLLogging, @SciMLMessage, AbstractVerbositySpecifier, AbstractVerbosityPreset,
    None, Minimal, Standard, Detailed, All, Silent, InfoLevel, WarnLevel, ErrorLevel,
    CustomLevel, AbstractMessageLevel, @verbosity_specifier

using SymbolicIndexingInterface: state_values, parameter_values

using ConcreteStructs: @concrete

import EnzymeCore

const CompiledFloats = Union{Float32, Float64}
import Preferences

abstract type AbstractNLSolverCache end
abstract type AbstractNLSolverAlgorithm end
abstract type AbstractNLSolver{algType, iip} end

function set_new_W! end
function set_W_γdt! end
function get_W end
function isfirstcall end
function isfirststage end
function isJcurrent end
function get_new_W_γdt_cutoff end
resize_J_W!(args...) = nothing
resize_nlsolver!(args...) = nothing

@enum MethodType begin
    DIRK
    COEFFICIENT_MULTISTEP
    NORDSIECK_MULTISTEP
    GLM
end

@enum NLStatus::Int8 begin
    FastConvergence = 2
    Convergence = 1
    SlowConvergence = 0
    VerySlowConvergence = -1
    Divergence = -2
end
const TryAgain = SlowConvergence

DEFAULT_PRECS(W, du, u, p, t, newW, Plprev, Prprev, solverdata) = nothing, nothing
isdiscretecache(cache) = false

@static if isdefined(DiffEqBase, :unitfulvalue)
    unitfulvalue(x) = DiffEqBase.unitfulvalue(x)
else
    unitfulvalue(x) = DiffEqBase.ForwardDiff.value(x)
end

include("doc_utils.jl")
include("misc_utils.jl")
include("verbosity.jl")

include("algorithms.jl")
include("composite_algs.jl")

include("alg_utils.jl")

include("caches/basic_caches.jl")

include("integrators/type.jl")
include("integrators/controllers.jl")
include("integrators/integrator_interface.jl")
include("integrators/integrator_utils.jl")
include("cache_utils.jl")
include("initialize_dae.jl")

include("perform_step/composite_perform_step.jl")

include("dense/generic_dense.jl")

include("iterator_interface.jl")
include("solve.jl")
include("initdt.jl")
include("interp_func.jl")

include("enzyme_rules.jl")

include("precompilation_setup.jl")

# Public API for shared loop functions and hooks.
# These are used by StochasticDiffEq to reuse ODE's time loop infrastructure
# while overriding SDE-specific behavior via hook dispatch.
public _savevalues!, _postamble!,
    handle_callbacks!, handle_tstop!, solution_endpoint_match_cur_integrator!,
    post_step_reject!, on_u_modified_at_init!, post_apply_step!,
    interp_at_saveat, post_savevalues!, finalize_solution_storage!,
    finalize_endpoint!, on_callbacks_complete!, is_composite_cache,
    is_composite_algorithm, final_progress

end
