module OrdinaryDiffEqNonlinearSolve

import ADTypes: AutoFiniteDiff, AutoForwardDiff

import SciMLBase
using SciMLBase: DAEFunction, DEIntegrator, NonlinearProblem, NonlinearLeastSquaresProblem, ODEProblem, DAEProblem
import DiffEqBase
import PreallocationTools
using SimpleNonlinearSolve: SimpleTrustRegion, SimpleGaussNewton
using NonlinearSolve: FastShortcutNonlinearPolyalg, FastShortcutNLLSPolyalg
using MuladdMacro, FastBroadcast
import FastClosures: @closure

import OrdinaryDiffEq: resize_nlsolver!, _initialize_dae!, OrdinaryDiffEqDifferentiation,
AbstractNLSolverAlgorithm, AbstractNLSolverCache, AbstractNLSolver, NewtonAlgorithm, @unpack,
OverrideInit, ShampineCollocationInit, BrownFullBasicInit, isnewton
import OrdinaryDiffEq.OrdinaryDiffEqDifferentiation: update_W!
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

include("type.jl")
include("utils.jl")
include("nlsolve.jl")
include("functional.jl")
include("newton.jl")
include("initialize_dae.jl")

end
