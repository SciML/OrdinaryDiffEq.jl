module OrdinaryDiffEqNonlinearSolve

using ADTypes: ADTypes, AutoForwardDiff, AutoFiniteDiff
using CommonSolve: init, solve, solve!, step!

import SciMLBase
import SciMLBase: remake
using SciMLOperators: update_coefficients!
using SciMLBase: DAEFunction, DEIntegrator, NonlinearFunction, NonlinearProblem,
    NonlinearLeastSquaresProblem, LinearProblem, ODEProblem, DAEProblem,
    get_tmp_cache, ReturnCode,
    AbstractNonlinearProblem, LinearAliasSpecifier,
    _vec, _reshape, postamble!, alg_order, isadaptive
import DiffEqBase
import DiffEqBase: OrdinaryDiffEqTag, calculate_residuals, calculate_residuals!,
    BrownFullBasicInit, ShampineCollocationInit, strip_domain_checks
import ConstructionBase
import PreallocationTools: DiffCache, get_tmp
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg, NewtonRaphson
# The operator Jacobian path is implemented in NonlinearSolveBase and needs its own floor.
import NonlinearSolveBase
# `get_u`/`get_fu` are the only inner-state reads that hold for every inner cache type:
# polyalgorithm caches keep `u`/`fu` on the active branch, not as top-level fields.
using NonlinearSolveBase:
    ArcLengthContinuation, HomotopyPolyAlgorithm, HomotopySweep, KantorovichHomotopy,
    get_linear_cache, get_u, get_fu
using MuladdMacro: @muladd
using FastBroadcast: @..
import FastClosures: @closure
using LinearAlgebra: UniformScaling, UpperTriangular, givens, cond, dot, lmul!, axpy!,
    I, rmul!, norm, mul!, ldiv!
import LinearAlgebra
import ArrayInterface: ArrayInterface
import LinearSolve
import ForwardDiff: ForwardDiff
using ForwardDiff: Dual
using RecursiveArrayTools: recursivecopy!
import OrdinaryDiffEqCore

import SciMLOperators: islinear, AbstractSciMLOperator, MatrixOperator
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
    alg_extrapolates, resize_J_W!, alg_autodiff,
    find_algebraic_vars_eqs

import OrdinaryDiffEqCore: _initialize_dae!,
    isnewton, get_W, isfirstcall, isfirststage,
    isJcurrent, get_new_W_γdt_cutoff, resize_nlsolver!, apply_step!,
    @SciMLMessage

import OrdinaryDiffEqDifferentiation: update_W!, is_always_new, build_uf, build_J_W,
    WOperator, StaticWOperator, wrapprecs, default_krylov_warm_start,
    build_jac_config, dolinsolve,
    resize_jac_config!, jacobian2W!, jacobian!

import StaticArraysCore: StaticArray

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("homotopy.jl")
include("initialize_dae.jl")

export BrownFullBasicInit, ShampineCollocationInit

# Nonlinear-solver algorithms accepted by OrdinaryDiffEq implicit methods and the
# solver-author hooks extended by sibling integrator packages.
# The `public` keyword is only parseable on Julia >= 1.11.0-DEV.469, so it is
# gated to keep the 1.10 floor parsing.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Expr(
            :public,
            :NLNewton, :NLFunctional, :NLAnderson, :HomotopyNonlinearSolveAlg,
            :NonlinearSolveAlg
        )
    )
end

end
