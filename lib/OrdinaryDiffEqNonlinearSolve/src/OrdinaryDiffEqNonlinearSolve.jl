module OrdinaryDiffEqNonlinearSolve

using ADTypes

import SciMLBase
import SciMLBase: init, solve, solve!, remake
using SciMLBase: DAEFunction, DEIntegrator, NonlinearFunction, NonlinearProblem,
                 NonlinearLeastSquaresProblem, LinearProblem, ODEProblem, DAEProblem,
                 update_coefficients!, get_tmp_cache, AbstractSciMLOperator, ReturnCode,
                 AbstractNonlinearProblem, LinearAliasSpecifier
import DiffEqBase
import PreallocationTools
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg, NewtonRaphson,
                      step!
using MuladdMacro, FastBroadcast
import FastClosures: @closure
using LinearAlgebra: UniformScaling, UpperTriangular, givens, cond, dot, lmul!, axpy!
import LinearAlgebra
import ArrayInterface
import LinearSolve
import ForwardDiff
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
                          AbstractNLSolver, NewtonAlgorithm, @unpack,
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
                           postamble!, isnonlinearsolve

import OrdinaryDiffEqDifferentiation: update_W!, is_always_new, build_uf, build_J_W,
                                      WOperator, StaticWOperator, wrapprecs,
                                      build_jac_config, dolinsolve, alg_autodiff,
                                      resize_jac_config!, do_newJW

import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

export BrownFullBasicInit, ShampineCollocationInit

end
