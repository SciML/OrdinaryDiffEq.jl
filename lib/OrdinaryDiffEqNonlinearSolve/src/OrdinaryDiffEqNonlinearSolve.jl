module OrdinaryDiffEqNonlinearSolve

using ADTypes: ADTypes, dense_ad, AutoForwardDiff, AutoFiniteDiff

import SciMLBase
import SciMLBase: init, solve, solve!, remake
using SciMLBase: DAEFunction, DEIntegrator, NonlinearFunction, NonlinearProblem,
    NonlinearLeastSquaresProblem, LinearProblem, ODEProblem, DAEProblem,
    update_coefficients!, get_tmp_cache, AbstractSciMLOperator, ReturnCode,
    AbstractNonlinearProblem, LinearAliasSpecifier
import DiffEqBase
import PreallocationTools: dualcache, get_tmp
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg, NewtonRaphson,
    step!, NonlinearVerbosity
using MuladdMacro: @muladd
using FastBroadcast: @..
import FastClosures: @closure
using LinearAlgebra: UniformScaling, UpperTriangular, givens, cond, dot, lmul!, axpy!
import LinearAlgebra
using SparseArrays: SparseMatrixCSC
import ArrayInterface: ArrayInterface, ismutable, restructure
import LinearSolve: OperatorAssumptions
import LinearSolve
import ForwardDiff: ForwardDiff, pickchunksize
using ForwardDiff: Dual
using LinearSolve: I, rmul!, norm, mul!, ldiv!
using RecursiveArrayTools: recursivecopy!
import SciMLStructures: canonicalize, Tunable, isscimlstructure
import OrdinaryDiffEqCore

import SciMLOperators: islinear
import OrdinaryDiffEqCore: nlsolve_f, set_new_W!, set_W_γdt!

@static if isdefined(OrdinaryDiffEqCore, :default_nlsolve)
    import OrdinaryDiffEqCore: default_nlsolve
end

using OrdinaryDiffEqCore: resize_nlsolver!, _initialize_dae!,
    AbstractNLSolverAlgorithm, AbstractNLSolverCache,
    AbstractNLSolver, NewtonAlgorithm,
    OverrideInit, ShampineCollocationInit, BrownFullBasicInit,
    _vec, _unwrap_val, DAEAlgorithm,
    _reshape, calculate_residuals, calculate_residuals!,
    has_special_newton_error, isadaptive,
    TryAgain, DIRK, COEFFICIENT_MULTISTEP, NORDSIECK_MULTISTEP, GLM,
    FastConvergence, Convergence,
    SlowConvergence, VerySlowConvergence, Divergence, NLStatus,
    MethodType, alg_order, error_constant,
    alg_extrapolates, resize_J_W!, has_autodiff

import OrdinaryDiffEqCore: _initialize_dae!, isnewton, get_W, isfirstcall, isfirststage,
    isJcurrent, get_new_W_γdt_cutoff, resize_nlsolver!, apply_step!,
    postamble!, @SciMLMessage

import OrdinaryDiffEqDifferentiation: update_W!, is_always_new, build_uf, build_J_W,
    WOperator, StaticWOperator, wrapprecs,
    build_jac_config, dolinsolve, alg_autodiff,
    resize_jac_config!,
    _JacConfigFiniteDiff

import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
    StaticMatrix

include("type.jl")
include("utils.jl")
include("deprecated.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

export BrownFullBasicInit, ShampineCollocationInit

end
