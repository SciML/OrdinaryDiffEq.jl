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

using LinearAlgebra: opnorm, I, UniformScaling, diag, isdiag

import PrecompileTools


import FastPower: fastpower

# Interfaces
import SciMLBase: solve!, step!, isadaptive
import DiffEqBase: initialize!

# DAE Initialization algorithms.
# `ShampineCollocationInit`/`BrownFullBasicInit` are unused here but re-exported
# for dependent OrdinaryDiffEq.jl sublibraries that import them from this package.
import DiffEqBase: DefaultInit, ShampineCollocationInit, BrownFullBasicInit

# Internal utils. `DEVerbosity` is re-exported for dependent sublibraries.
import DiffEqBase: ODE_DEFAULT_NORM,
    ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE,
    ODE_DEFAULT_UNSTABLE_CHECK,
    DEVerbosity, _process_verbose_param

import SciMLOperators: MatrixOperator, FunctionOperator,
    update_coefficients, update_coefficients!,
    isconstant

import Random

import RecursiveArrayTools: recursivecopy!, recursivecopy, recursive_bottom_eltype, recursive_unitless_bottom_eltype, recursive_unitless_eltype, copyat_or_push!, DiffEqArray

import RecursiveArrayTools
using BinaryHeaps: BinaryHeap, FasterForward
using ArrayInterface: ArrayInterface

import TruncatedStacktraces: @truncate_stacktrace, VERBOSE_MSG

# Integrator Interface
import Base: resize!
import SciMLBase: deleteat!, addat!, full_cache, user_cache, u_cache, du_cache,
    resize_non_user_cache!, deleteat_non_user_cache!, addat_non_user_cache!,
    terminate!, get_du, get_dt, get_proposed_dt, set_proposed_dt!,
    savevalues!,
    add_tstop!, has_tstop, first_tstop, pop_tstop!,
    add_saveat!, set_reltol!,
    set_abstol!, postamble!, last_step_failed
import DiffEqBase: get_tstops, get_tstops_array

# `check_error!` is owned by (and public in) SciMLBase, re-exported through DiffEqBase.
# `_vec`/`_reshape`/`unwrap_cache` are owned by SciMLBase and re-exported here for
# dependent OrdinaryDiffEq.jl sublibraries that import them from this package.
using SciMLBase: check_error!, _vec, _reshape, unwrap_cache

using FastBroadcast: @..

using FunctionWrappers: FunctionWrapper

using SciMLBase: NoInit, CheckInit, OverrideInit, AbstractDEProblem, _unwrap_val,
    ODEAliasSpecifier

# Names previously relied on implicitly through `@reexport using SciMLBase`.
using SciMLBase: SciMLBase, CallbackSet, ContinuousCallback, DAEProblem,
    DAESolution, DiscreteProblem, DynamicalODEFunction,
    IntervalNonlinearProblem, NonlinearLeastSquaresProblem,
    ODEProblem, ReturnCode, SplitFunction, SplitSDEFunction,
    VectorContinuousCallback, auto_dt_reset!, derivative_discontinuity!,
    get_tmp_cache, isinplace, reinit!
using SciMLOperators: SciMLOperators
using CommonSolve: solve

import SciMLBase: AbstractNonlinearProblem, alg_order, LinearAliasSpecifier

import SciMLBase: islinear
# `calculate_residuals`/`calculate_residuals!` are unused here but re-exported for
# dependent OrdinaryDiffEq.jl sublibraries that import them from this package.
import DiffEqBase: timedepentdtmin, calculate_residuals, calculate_residuals!

# MacroTools and Adapt imported but not directly used in OrdinaryDiffEqCore
# using MacroTools, Adapt
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AutoSparse, dense_ad
import Accessors: @reset

# SciMLStructures symbols imported but not directly used in OrdinaryDiffEqCore
# using SciMLStructures: canonicalize, Tunable, isscimlstructure

# `Minimal` is unused here but re-exported for dependent OrdinaryDiffEq.jl
# sublibraries that import it from this package.
using SciMLLogging: SciMLLogging, @SciMLMessage, Standard, Minimal

using SymbolicIndexingInterface: state_values
import SymbolicIndexingInterface: parameter_values

using EnumX: @enumx

import EnzymeCore

# Per-stage Newton initial-guess ("predictor") strategies for implicit RK methods.
@enumx Predictor begin
    Trivial        # zero increment (z = 0)
    Linear         # linear extrapolation (z = dt * fsalfirst)
    MaxOrder       # full previous-step interpolant
    VariableOrder  # interpolant, order reduced as the stage extrapolates further
    CutoffOrder    # interpolant, full order below a cutoff abscissa else order 1
    CopyPrev       # reuse the previous stage's derivative
    StageExtrap    # extrapolate recent stage derivatives
    Tableau        # tableau-derived predictor (α / const_stage_guess)
end
export Predictor

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

isdiscretecache(cache) = false

unitfulvalue(x) = SciMLBase.unitfulvalue(x)

# Declare the documented custom-stepsize-controller author interface `public` so downstream
# packages can drop their `OrdinaryDiffEqCore.X` non-public ExplicitImports ignores. These
# are the documented "Required methods" of the controller API (docs/src/api/controllers.md).
# The `public` keyword only parses on Julia >= 1.11.0-DEV.469, so it is gated to keep the
# 1.10 floor parsing.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Expr(
            :public,
            :stepsize_controller!, :step_accept_controller!, :step_reject_controller!
        )
    )
end

include("doc_utils.jl")
include("misc_utils.jl")

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
include("disco.jl")

include("dense/generic_dense.jl")

include("iterator_interface.jl")
include("solve.jl")
include("initdt.jl")
include("interp_func.jl")

include("enzyme_rules.jl")

include("precompilation_setup.jl")

end
