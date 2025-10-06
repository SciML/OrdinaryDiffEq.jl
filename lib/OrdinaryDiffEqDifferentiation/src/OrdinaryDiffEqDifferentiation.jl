module OrdinaryDiffEqDifferentiation

import ADTypes
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AbstractADType, AutoSparse

import ForwardDiff, FiniteDiff
import ForwardDiff.Dual
import LinearSolve
import LinearSolve: OperatorAssumptions
import FunctionWrappersWrappers
import DiffEqBase

import LinearAlgebra
import LinearAlgebra: Diagonal, I, UniformScaling, diagind, mul!, lmul!, axpby!, opnorm, lu
import LinearAlgebra: LowerTriangular, UpperTriangular
import SparseArrays: SparseMatrixCSC, AbstractSparseMatrix, nonzeros, sparse, spzeros
import ArrayInterface
import ArrayInterface: fast_scalar_indexing, zeromatrix, lu_instance

# StaticArrayInterface imported but not used
# import StaticArrayInterface
import StaticArrays
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

using DiffEqBase: TimeGradientWrapper,
                  UJacobianWrapper, TimeDerivativeWrapper,
                  UDerivativeWrapper
import SciMLBase: SciMLBase, constructorof, @set, isinplace, has_jvp, unwrapped_f, DEIntegrator, ODEFunction, SplitFunction, DynamicalODEFunction, DAEFunction, islinear, remake, solve!, isconstant
using SciMLBase: @set, @reset
import SciMLOperators: SciMLOperators, IdentityOperator, update_coefficients, update_coefficients!, MatrixOperator, AbstractSciMLOperator, ScalarOperator
import SparseMatrixColorings: ConstantColoringAlgorithm, GreedyColoringAlgorithm, ColoringProblem,
                               ncolors, column_colors, coloring, sparsity_pattern
import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                          DAEAlgorithm,
                          OrdinaryDiffEqImplicitAlgorithm, CompositeAlgorithm,
                          OrdinaryDiffEqExponentialAlgorithm,
                          OrdinaryDiffEqAdaptiveExponentialAlgorithm, @unpack,
                          AbstractNLSolver, nlsolve_f, issplit,
                          concrete_jac, unwrap_alg, OrdinaryDiffEqCache, _vec, standardtag,
                          isnewton, _unwrap_val,
                          set_new_W!, set_W_γdt!, alg_difftype, unwrap_cache, diffdir,
                          get_W, isfirstcall, isfirststage, isJcurrent,
                          get_new_W_γdt_cutoff,
                          TryAgain, DIRK, COEFFICIENT_MULTISTEP, NORDSIECK_MULTISTEP, GLM,
                          FastConvergence, Convergence, SlowConvergence,
                          VerySlowConvergence, Divergence, NLStatus, MethodType, constvalue

import OrdinaryDiffEqCore: get_chunksize, resize_J_W!, resize_nlsolver!, alg_autodiff,
                           _get_fwd_tag

import ConstructionBase
using ConstructionBase: constructorof

import DifferentiationInterface as DI

using FastBroadcast: @..

using ConcreteStructs: @concrete

@static if isdefined(SciMLBase, :OrdinaryDiffEqTag)
    import SciMLBase: OrdinaryDiffEqTag
elseif isdefined(DiffEqBase, :OrdinaryDiffEqTag)
    import DiffEqBase: OrdinaryDiffEqTag
else
    struct OrdinaryDiffEqTag end
end

include("alg_utils.jl")
include("linsolve_utils.jl")
include("derivative_utils.jl")
include("derivative_wrappers.jl")
include("operators.jl")

end
