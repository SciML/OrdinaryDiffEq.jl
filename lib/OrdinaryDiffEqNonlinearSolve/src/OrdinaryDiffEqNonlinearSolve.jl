module OrdinaryDiffEqNonlinearSolve

import ADTypes: AutoFiniteDiff, AutoForwardDiff

import SciMLBase
import SciMLBase: init, solve, solve!
using SciMLBase: DAEFunction, DEIntegrator, NonlinearProblem, NonlinearLeastSquaresProblem, LinearProblem, ODEProblem, DAEProblem, update_coefficients!, get_tmp_cache, AbstractSciMLOperator
import DiffEqBase
import PreallocationTools
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg
using MuladdMacro, FastBroadcast
import FastClosures: @closure
using LinearAlgebra: UniformScaling
import ArrayInterface
import LinearSolve
using LinearSolve: I

import SciMLOperators: islinear
import OrdinaryDiffEq: nlsolve_f, set_new_W!, set_W_γdt!

using OrdinaryDiffEq: OrdinaryDiffEqDifferentiation,
AbstractNLSolverAlgorithm, AbstractNLSolverCache, AbstractNLSolver, NewtonAlgorithm, @unpack,
OverrideInit, ShampineCollocationInit, BrownFullBasicInit, _vec, _unwrap_val, DAEAlgorithm,
_reshape, calculate_residuals, calculate_residuals!, has_special_newton_error, isadaptive

import OrdinaryDiffEq: _initialize_dae!, resize_nlsolver!, isnewton, get_W, isfirstcall

import OrdinaryDiffEq.OrdinaryDiffEqDifferentiation: update_W!, is_always_new, build_uf, build_J_W, WOperator, StaticWOperator, wrapprecs, build_jac_config, dolinsolve
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

end
