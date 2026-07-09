module OrdinaryDiffEqNonlinearSolve

using ADTypes: ADTypes, AutoForwardDiff, AutoFiniteDiff

import SciMLBase
import SciMLBase: init, solve, remake
using SciMLBase: DAEFunction, DEIntegrator, NonlinearFunction, NonlinearProblem,
    NonlinearLeastSquaresProblem, LinearProblem, ODEProblem, DAEProblem,
    update_coefficients!, get_tmp_cache, ReturnCode,
    AbstractNonlinearProblem, LinearAliasSpecifier,
    _vec, _reshape, postamble!, alg_order, isadaptive
import DiffEqBase
import DiffEqBase: OrdinaryDiffEqTag, calculate_residuals, calculate_residuals!,
    BrownFullBasicInit, ShampineCollocationInit
import ConstructionBase
import PreallocationTools: DiffCache, get_tmp
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg, NewtonRaphson,
    step!
using MuladdMacro: @muladd
using FastBroadcast: @..
import FastClosures: @closure
using LinearAlgebra: UniformScaling, UpperTriangular, givens, cond, dot, lmul!, axpy!,
    I, rmul!, norm, mul!, ldiv!
import LinearAlgebra
using SparseArrays: SparseMatrixCSC
import ArrayInterface: ArrayInterface
import LinearSolve
import ForwardDiff: ForwardDiff
using ForwardDiff: Dual
using RecursiveArrayTools: recursivecopy!
import OrdinaryDiffEqCore

import SciMLOperators: islinear, AbstractSciMLOperator
import OrdinaryDiffEqCore: nlsolve_f, set_new_W!, set_W_γdt!

import OrdinaryDiffEqCore: default_nlsolve

using OrdinaryDiffEqCore: resize_nlsolver!, _initialize_dae!,
    AbstractNLSolverAlgorithm, AbstractNLSolverCache,
    AbstractNLSolver, NewtonAlgorithm,
    DAEAlgorithm,
    has_special_newton_error,
    TryAgain, DIRK, COEFFICIENT_MULTISTEP,
    Convergence,
    Divergence, NLStatus,
    MethodType, error_constant,
    alg_extrapolates, resize_J_W!, alg_autodiff

import OrdinaryDiffEqCore: _initialize_dae!,
    isnewton, get_W, isfirstcall, isfirststage,
    isJcurrent, get_new_W_γdt_cutoff, resize_nlsolver!, apply_step!,
    @SciMLMessage

import OrdinaryDiffEqDifferentiation: update_W!, is_always_new, build_uf, build_J_W,
    WOperator, StaticWOperator, wrapprecs,
    build_jac_config, dolinsolve,
    resize_jac_config!, jacobian2W!, jacobian!

import StaticArraysCore: StaticArray

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

export BrownFullBasicInit, ShampineCollocationInit

# Declare the documented solver-author extension API `public` so downstream solver
# packages (the OrdinaryDiffEq solver sublibs, DelayDiffEq's fixed-point solver, etc.)
# can drop their `OrdinaryDiffEqNonlinearSolve.X` non-public ExplicitImports ignores.
# These are the nonlinear-solver algorithm types and the stepping/build hooks that
# downstream packages legitimately construct and extend. The `public` keyword is only
# parseable on Julia >= 1.11.0-DEV.469, so it is gated to keep the 1.10 floor parsing.
# Pure internal codegen/perf helpers are intentionally left non-public.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Expr(
            :public,
            # Nonlinear-solver algorithm types
            :NLNewton, :NLFunctional, :NLAnderson,
            # Solver construction
            :build_nlsolver,
            # Stepping / status hooks extended and called by downstream solvers
            :nlsolve!, :nlsolvefail, :compute_step!, :initial_η,
            :markfirststage!, :du_alias_or_new,
            # Anderson-acceleration kernels (DelayDiffEq fixed-point solver reuses these)
            :anderson, :anderson!
        )
    )
end

end
