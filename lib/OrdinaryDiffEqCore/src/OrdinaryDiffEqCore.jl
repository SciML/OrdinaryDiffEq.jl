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

using LinearAlgebra: opnorm, I, UniformScaling, diag, isdiag, Diagonal

import PrecompileTools


import FastPower: fastpower

import FindFirstFunctions: StrategyKind, KIND_BRACKET_GALLOP, KIND_INTERPOLATION_SEARCH,
    searchsorted_first, searchsorted_last

# Interfaces
import SciMLBase: solve!, step!, isadaptive
import DiffEqBase: initialize!

# DAE Initialization algorithms.
# `ShampineCollocationInit`/`BrownFullBasicInit` are unused here but re-exported
# for dependent OrdinaryDiffEq.jl sublibraries that import them from this package.
import DiffEqBase: DefaultInit, ShampineCollocationInit, BrownFullBasicInit

# Specialization level owned by SciMLBase (declared `public` there, not
# exported). Re-exported here so `using OrdinaryDiffEq` surfaces it unqualified.
import SciMLBase: AutoDePSpecialize
export AutoDePSpecialize

# Internal utils. `DEVerbosity` is re-exported for dependent sublibraries.
import DiffEqBase: ODE_DEFAULT_NORM,
    ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE,
    ODE_DEFAULT_UNSTABLE_CHECK,
    DEVerbosity, _process_verbose_param

import SciMLOperators: MatrixOperator, FunctionOperator,
    update_coefficients, update_coefficients!,
    isconstant

import Random
import Printf: @sprintf

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

import SciMLBase: AbstractNonlinearProblem, alg_order, LinearAliasSpecifier, log_numerical_instability, has_mtk_sys

import SciMLBase: islinear
# `calculate_residuals`/`calculate_residuals!` are unused here but re-exported for
# dependent OrdinaryDiffEq.jl sublibraries that import them from this package.
import DiffEqBase: timedepentdtmin, calculate_residuals, calculate_residuals!

# MacroTools and Adapt imported but not directly used in OrdinaryDiffEqCore
# using MacroTools, Adapt
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AutoSparse, dense_ad
import Accessors: @reset
import ConstructionBase

# SciMLStructures symbols imported but not directly used in OrdinaryDiffEqCore
# using SciMLStructures: canonicalize, Tunable, isscimlstructure

# `Minimal` is unused here but re-exported for dependent OrdinaryDiffEq.jl
# sublibraries that import it from this package.
using SciMLLogging: SciMLLogging, @SciMLMessage, Standard, Minimal

using SymbolicIndexingInterface: state_values
import SymbolicIndexingInterface: parameter_values

using EnumX: @enumx

import EnzymeCore

"""
    Predictor

Enumeration of per-stage initial-guess strategies for implicit Runge-Kutta
methods.

# Values
- `Predictor.Trivial`: Use a zero increment.
- `Predictor.Linear`: Use a linear extrapolation from the first-same-as-last
  derivative.
- `Predictor.MaxOrder`: Use the full previous-step interpolation order.
- `Predictor.VariableOrder`: Reduce interpolation order as the stage extrapolates
  farther from the previous step.
- `Predictor.CutoffOrder`: Use full interpolation order below a cutoff and order
  one above it.
- `Predictor.CopyPrev`: Reuse the previous stage derivative.
- `Predictor.StageExtrap`: Extrapolate recent stage derivatives.
- `Predictor.Tableau`: Use the tableau-derived stage guess.
"""
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

"""
    CompiledFloats

`Union{Float32, Float64}` — the floating-point element types for which the solvers
provide fully precompiled specializations.
"""
const CompiledFloats = Union{Float32, Float64}
import Preferences

"""
    AbstractNLSolverCache

Abstract supertype of the mutable cache held by an [`AbstractNLSolver`](@ref),
storing the W-matrix, factorizations, residual buffers, and per-iteration state.
"""
abstract type AbstractNLSolverCache end
"""
    AbstractNLSolverAlgorithm

Abstract supertype of the algorithm objects that configure a nonlinear solver
(e.g. `NLNewton`, `NLFunctional`, `NLAnderson`). Passed as the `nlsolve` keyword
of an implicit algorithm and consumed internally by `OrdinaryDiffEqNonlinearSolve.build_nlsolver`.
"""
abstract type AbstractNLSolverAlgorithm end
"""
    AbstractNLSolver{algType, iip}

Abstract supertype of the nonlinear solver object that implicit algorithms use to
solve their implicit stage equations. Concrete subtypes (e.g. `NLSolver` in
OrdinaryDiffEqNonlinearSolve) are built with `OrdinaryDiffEqNonlinearSolve.build_nlsolver`
and driven with `OrdinaryDiffEqNonlinearSolve.nlsolve!`. The type parameters are
the nonlinear-solver algorithm type and the in-place flag `iip`.
"""
abstract type AbstractNLSolver{algType, iip} end

"""
    set_new_W!(nlsolver, val::Bool) -> Bool

Set the flag recording whether a fresh `W` was just computed for this step and
return `val`. Read via `get_new_W!` to decide whether to reset the Newton
convergence estimate.
"""
function set_new_W! end
"""
    set_W_γdt!(nlsolver, W_γdt)

Store the `γΔt` value at which the current `W` was formed and return it. A change
in `γΔt` beyond [`get_new_W_γdt_cutoff`](@ref) triggers a `W` refactorization.
"""
function set_W_γdt! end
"""
    get_W(nlsolver)

Return the `W = M/(γΔt) - J` matrix (or its factorization) held by the nonlinear
solver's cache.
"""
function get_W end
"""
    isfirstcall(nlsolver) -> Bool

Return whether this is the first nonlinear solve of the current step (so a fresh
`W`/initial guess is needed).
"""
function isfirstcall end
"""
    isfirststage(nlsolver) -> Bool

Return whether the nonlinear solver is on the first implicit stage of a multi-stage
method (used to decide predictor/W reuse).
"""
function isfirststage end
"""
    isJcurrent(nlsolver, integrator) -> Bool

Return whether the Jacobian stored on the solver is current for the present
`integrator` state (so it need not be recomputed).
"""
function isJcurrent end
"""
    get_new_W_γdt_cutoff(nlsolver)

Return the relative-change threshold on `γΔt` above which the nonlinear solver
recomputes/refactorizes `W` rather than reusing the existing one.
"""
function get_new_W_γdt_cutoff end
"""
    resize_J_W!(nlsolver, integrator, i)

Resize the Jacobian and `W` matrices held by `nlsolver` to length `i` after the
state size changes (e.g. from a `resize!` callback). No-op fallback.
"""
resize_J_W!(args...) = nothing
"""
    resize_nlsolver!(integrator, i)

Resize the nonlinear solver's internal buffers to state length `i` after the state
size changes. No-op fallback.
"""
resize_nlsolver!(args...) = nothing

"""
    MethodType

`@enum` classifying how an implicit algorithm forms its `W = M/(γΔt) - J` matrix
and stage system. One of [`DIRK`](@ref), [`COEFFICIENT_MULTISTEP`](@ref),
[`NORDSIECK_MULTISTEP`](@ref), or [`GLM`](@ref). The nonlinear solver uses it to
scale `γW` appropriately in `OrdinaryDiffEqNonlinearSolve.nlsolve!`.
"""
@enum MethodType begin
    DIRK
    COEFFICIENT_MULTISTEP
    NORDSIECK_MULTISTEP
    GLM
end

"""
    DIRK

[`MethodType`](@ref) value for diagonally-implicit Runge–Kutta stage systems.
"""
DIRK

"""
    COEFFICIENT_MULTISTEP

[`MethodType`](@ref) value for coefficient-form multistep (e.g. BDF) methods.
"""
COEFFICIENT_MULTISTEP

"""
    NORDSIECK_MULTISTEP

[`MethodType`](@ref) value for Nordsieck-form multistep methods.
"""
NORDSIECK_MULTISTEP

"""
    GLM

[`MethodType`](@ref) value for general linear methods.
"""
GLM

"""
    NLStatus

`@enum` reporting the outcome/convergence quality of a nonlinear solve, ordered
from best to worst: [`FastConvergence`](@ref) (`2`), [`Convergence`](@ref) (`1`),
[`SlowConvergence`](@ref) (`0`), [`VerySlowConvergence`](@ref) (`-1`),
[`Divergence`](@ref) (`-2`). A non-positive value means the solve failed
(`OrdinaryDiffEqNonlinearSolve.nlsolvefail`).
"""
@enum NLStatus::Int8 begin
    FastConvergence = 2
    Convergence = 1
    SlowConvergence = 0
    VerySlowConvergence = -1
    Divergence = -2
end

"""
    FastConvergence

[`NLStatus`](@ref) value (`2`) indicating the nonlinear solve converged quickly.
"""
FastConvergence

"""
    Convergence

[`NLStatus`](@ref) value (`1`) indicating the nonlinear solve converged.
"""
Convergence

"""
    SlowConvergence

[`NLStatus`](@ref) value (`0`) indicating slow convergence. Also aliased as
[`TryAgain`](@ref).
"""
SlowConvergence

"""
    VerySlowConvergence

[`NLStatus`](@ref) value (`-1`) indicating very slow convergence; treated as a
failure.
"""
VerySlowConvergence

"""
    Divergence

[`NLStatus`](@ref) value (`-2`) indicating the nonlinear iteration diverged.
"""
Divergence

"""
    TryAgain

Alias for [`SlowConvergence`](@ref); a sentinel [`NLStatus`](@ref) signalling that
the step should be retried (e.g. with a fresh Jacobian).
"""
const TryAgain = SlowConvergence

"""
    isdiscretecache(cache) -> Bool

Return whether `cache` belongs to a discrete-time / map-iteration algorithm rather
than a continuous ODE solver (`false` by default). Companion to
[`isdiscretealg`](@ref).
"""
isdiscretecache(cache) = false

"""
    unitfulvalue(x)

Return the unit-stripped numeric value of `x` (delegates to
`SciMLBase.unitfulvalue`). Used where a bare number is needed from a possibly
`Unitful` quantity.
"""
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

# Solver-author extension API. These names are the documented/intended interface
# that downstream OrdinaryDiffEq.jl / StochasticDiffEq.jl solver packages subtype,
# extend, or call. They are made public (not exported) so that ExplicitImports'
# public-API checks recognize them as the supported extension surface. Genuine
# codegen/perf internals (@fold/@threaded/@OnDemandTableauExtract/@swap!) and
# precompile-workload helpers are deliberately NOT included here.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Expr(
            :public,
            :AbstractController, :AbstractControllerCache, :AbstractNLSolver, :AbstractNLSolverAlgorithm, :AbstractThreadingOption,
            :accept_step_controller, :alg_adaptive_order, :alg_autodiff, :alg_cache, :alg_difftype, :alg_extrapolates,
            :alg_maximum_order, :alg_stability_size, :AutoAlgSwitch, :AutoSwitch, :BaseThreads,
            :beta1_default, :beta2_default, Symbol("@cache"), :COEFFICIENT_MULTISTEP, :CommonControllerOptions, :CompositeAlgorithm,
            :CompositeController, :constvalue, :Convergence, :current_extrapolant,
            :current_interpolant, :DAEAlgorithm, :default_autoswitch, :default_controller, :default_linear_interpolation,
            :default_nlsolve, :DEOptions, :DIRK, :Divergence, :dt_required, :DummyController,
            :explicit_rk_docstring, :ExponentialAlgorithm, :FastConvergence, :gamma_default, :generic_solver_docstring,
            :get_current_adaptive_order, :get_current_alg_autodiff, :get_differential_vars, :get_EEst, :get_failfactor, :get_fsalfirstlast,
            :get_gamma, :get_new_W_γdt_cutoff, :get_qmax, :get_qmax_first_step, :get_qmin, :get_qsteady_max,
            :get_qsteady_min, :get_W, :GLM, :hermite_interpolant, :IController,
            :ImplicitSecondOrderAlgorithm, :increment_accept!, :increment_nf!, :increment_reject!, :InterpolationData, :isautoswitch,
            :is_composite_algorithm, :isdefaultalg, :isdiscretecache, :isdtchangeable, :isfirstcall,
            :isfirststage, :isfsal, :isimplicit, :isJcurrent, :is_mass_matrix_alg, :ismultistep,
            :issplit, :isthreaded, :isWmethod, :MethodType, :NewtonAlgorithm, :nlsolve_f,
            :NLStatus, :NORDSIECK_MULTISTEP, :_ode_addsteps!, :ode_addsteps!, :ODEIntegrator, :_ode_interpolant,
            :OrdinaryDiffEqAdamsVarOrderVarStepAlgorithm, :OrdinaryDiffEqAdaptiveAlgorithm, :OrdinaryDiffEqAdaptiveExponentialAlgorithm, :OrdinaryDiffEqAdaptiveImplicitAlgorithm, :OrdinaryDiffEqAdaptiveImplicitSecondOrderAlgorithm, :OrdinaryDiffEqAdaptivePartitionedAlgorithm,
            :OrdinaryDiffEqAlgorithm, :OrdinaryDiffEqCache, :OrdinaryDiffEqCompositeAlgorithm, :OrdinaryDiffEqConstantCache, :OrdinaryDiffEqExponentialAlgorithm, :OrdinaryDiffEqImplicitAlgorithm,
            :OrdinaryDiffEqImplicitSecondOrderAlgorithm, :OrdinaryDiffEqInterpolation, :OrdinaryDiffEqLinearExponentialAlgorithm, :OrdinaryDiffEqMutableCache, :OrdinaryDiffEqNewtonAdaptiveAlgorithm, :OrdinaryDiffEqNewtonAlgorithm,
            :OrdinaryDiffEqPartitionedAlgorithm, :OrdinaryDiffEqRosenbrockAdaptiveAlgorithm, :OrdinaryDiffEqRosenbrockAlgorithm, :PartitionedAlgorithm, :perform_step!, :PIController,
            :PIDController, :PolyesterThreads, :post_newton_controller!, :PredictiveController,
            :qmax_default, :qmin_default, :reinit_controller!, :resize_J_W!, :resize_nlsolver!,
            :RosenbrockAlgorithm, :Sequential, :set_EEst!, :set_new_W!, :setup_controller_cache, :set_W_γdt!,
            :SlowConvergence, :step_accept_controller!, :step_reject_controller!, :stepsize_controller!, :StochasticDiffEqAdaptiveAlgorithm, :StochasticDiffEqAlgorithm,
            :StochasticDiffEqCompositeAlgorithm, :StochasticDiffEqJumpAdaptiveAlgorithm, :StochasticDiffEqJumpAlgorithm, :StochasticDiffEqJumpDiffusionAdaptiveAlgorithm,
            :StochasticDiffEqJumpDiffusionAlgorithm, :StochasticDiffEqJumpNewtonAdaptiveAlgorithm, :StochasticDiffEqJumpNewtonDiffusionAdaptiveAlgorithm, :StochasticDiffEqNewtonAdaptiveAlgorithm, :StochasticDiffEqNewtonAlgorithm,
            :StochasticDiffEqRODEAdaptiveAlgorithm, :StochasticDiffEqRODEAlgorithm, :StochasticDiffEqRODECompositeAlgorithm, :sync_controllers!, :TryAgain,
            :unwrap_alg, :uses_uprev, :VerySlowConvergence,
            # Round 2: remaining cross-sublib / extension surface owned by OrdinaryDiffEqCore.
            # Error / sentinel / dispatch-helper types shared across solver sublibs.
            :CompiledFloats, :DerivativeOrderNotPossibleError, :DifferentialVarsUndefined,
            # Interpolation kernels (companions to the public _ode_interpolant / ode_addsteps!).
            :_ode_interpolant!, :ode_interpolant, :ode_interpolant!, :hermite_interpolant!,
            :current_extrapolant!, :interpolation_differential_vars,
            # Algorithm-trait predicates extended/queried by solver sublibs.
            :has_stage_limiter,
            :standardtag, :concrete_jac, :has_autodiff, :has_dtnew_modification,
            :has_special_newton_error, :has_stiff_interpolation, :alg_can_repeat_jac,
            :allows_null_u0, :isaposteriori, :isdiscretealg, :isdp8,
            :isesdirk, :isfirk, :isnewton, :only_diagonal_mass_matrix, :fsal_typeof,
            :ssp_coefficient, :fac_default_gamma, :qsteady_max_default, :qsteady_min_default,
            # Order / stepsize / autodiff-config accessors used across sublibs.
            :get_current_alg_order, :get_current_qmax, :get_chunksize, :_get_fdtype,
            :_get_fwd_chunksize, :_get_fwd_chunksize_int, :_fixup_ad, :diffdir,
            :error_constant, :unitfulvalue,
            # Integrator step / cache / initialization hooks.
            :_ode_init, :_determine_initdt, :ode_determine_initdt, :_initialize_dae!,
            :find_algebraic_vars_eqs, :postamble!, :apply_step!, :last_step_failed, :reset_alg_dependent_opts!,
            :set_discontinuity, :resolve_basic,
            # Noise hooks used by the SDE/RODE solver sublibs.
            :accept_noise!, :reinit_noise!, :reject_noise!, :save_noise!, :noise_curt,
            :is_noise_saveable,
            # Docstring builder used by solver sublibs.
            :differentiation_rk_docstring,
        )
    )
end

end
