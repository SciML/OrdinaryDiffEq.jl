module OrdinaryDiffEqNonlinearSolve

import ADTypes: AutoFiniteDiff, AutoForwardDiff

import SciMLBase
using SciMLBase: DAEFunction, DEIntegrator, NonlinearProblem, NonlinearLeastSquaresProblem, ODEProblem, DAEProblem, update_coefficients!, get_tmp_cache, AbstractSciMLOperator
import DiffEqBase
import PreallocationTools
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg
using MuladdMacro, FastBroadcast
import FastClosures: @closure
using LinearAlgebra: UniformScaling
import ArrayInterface
using LinearSolve: I

import OrdinaryDiffEq: nlsolve_f, set_new_W!, set_W_γdt!

using OrdinaryDiffEq: resize_nlsolver!, _initialize_dae!, OrdinaryDiffEqDifferentiation,
AbstractNLSolverAlgorithm, AbstractNLSolverCache, AbstractNLSolver, NewtonAlgorithm, @unpack,
OverrideInit, ShampineCollocationInit, BrownFullBasicInit, _vec, _unwrap_val, DAEAlgorithm,
_reshape, calculate_residuals, calculate_residuals!, has_special_newton_error, isadaptive

import OrdinaryDiffEq: isnewton

import OrdinaryDiffEq.OrdinaryDiffEqDifferentiation: update_W!, is_always_new, build_uf, build_J_W, WOperator, StaticWOperator
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

end
