module OrdinaryDiffEqCore

if isdefined(Base, :Experimental) &&
   isdefined(Base.Experimental, Symbol("@max_methods"))
    @eval Base.Experimental.@max_methods 1
end

using DocStringExtensions
using Reexport
@reexport using DiffEqBase

using Logging

using MuladdMacro, FastClosures

using LinearAlgebra

using PrecompileTools

import FillArrays: Trues, Falses

import FastPower

# Interfaces
import DiffEqBase: solve!, step!, initialize!, isadaptive

# Internal utils
import DiffEqBase: ODE_DEFAULT_NORM,
                   ODE_DEFAULT_ISOUTOFDOMAIN, ODE_DEFAULT_PROG_MESSAGE,
                   ODE_DEFAULT_UNSTABLE_CHECK

import SciMLOperators: AbstractSciMLOperator, AbstractSciMLScalarOperator,
                       MatrixOperator, FunctionOperator,
                       update_coefficients, update_coefficients!, DEFAULT_UPDATE_FUNC,
                       isconstant

using DiffEqBase: DEIntegrator

import RecursiveArrayTools: chain, recursivecopy!

using SimpleUnPack, RecursiveArrayTools, DataStructures, ArrayInterface

import TruncatedStacktraces

import StaticArraysCore: SArray, MVector, SVector, StaticArray, MMatrix,
                         StaticMatrix

# Integrator Interface
import DiffEqBase: resize!, deleteat!, addat!, full_cache, user_cache, u_cache, du_cache,
                   resize_non_user_cache!, deleteat_non_user_cache!, addat_non_user_cache!,
                   terminate!, get_du, get_dt, get_proposed_dt, set_proposed_dt!,
                   u_modified!, savevalues!,
                   add_tstop!, has_tstop, first_tstop, pop_tstop!,
                   add_saveat!, set_reltol!,
                   set_abstol!, postamble!, last_step_failed,
                   isautodifferentiable,
                   get_tstops, get_tstops_array, get_tstops_max

using DiffEqBase: check_error!, @def, _vec, _reshape

using FastBroadcast: @.., True, False

using SciMLBase: NoInit, CheckInit, OverrideInit, AbstractDEProblem, _unwrap_val,
                 ODEAliasSpecifier

import SciMLBase: AbstractNonlinearProblem, alg_order, LinearAliasSpecifier

import DiffEqBase: calculate_residuals,
                   calculate_residuals!, unwrap_cache,
                   @tight_loop_macros,
                   islinear, timedepentdtmin

import Polyester
using MacroTools, Adapt
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AbstractADType
import Accessors: @reset

using SciMLStructures: canonicalize, Tunable, isscimlstructure

using SymbolicIndexingInterface: state_values, parameter_values, is_variable,
                                 variable_index,
                                 symbolic_type, NotSymbolic

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

include("precompilation_setup.jl")

end
