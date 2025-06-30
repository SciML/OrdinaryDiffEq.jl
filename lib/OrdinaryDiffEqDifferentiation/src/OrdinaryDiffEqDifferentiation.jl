module OrdinaryDiffEqDifferentiation

import ADTypes
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AbstractADType, AutoSparse

import ForwardDiff, FiniteDiff
import ForwardDiff.Dual
import LinearSolve
import LinearSolve: OperatorAssumptions
import FunctionWrappersWrappers
using DiffEqBase

import LinearAlgebra
import LinearAlgebra: Diagonal, I, UniformScaling, diagind, mul!, lmul!, axpby!, opnorm, lu
import LinearAlgebra: LowerTriangular, UpperTriangular
import SparseArrays: SparseMatrixCSC, AbstractSparseMatrix, nonzeros, sparse, spzeros
import ArrayInterface

import StaticArrayInterface
import StaticArrays
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

using DiffEqBase: TimeGradientWrapper,
                  UJacobianWrapper, TimeDerivativeWrapper,
                  UDerivativeWrapper
using SciMLBase: AbstractSciMLOperator, constructorof, @set
using SciMLOperators
import SparseMatrixColorings
import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm,
                          DAEAlgorithm,
                          OrdinaryDiffEqImplicitAlgorithm, CompositeAlgorithm,
                          OrdinaryDiffEqExponentialAlgorithm,
                          OrdinaryDiffEqAdaptiveExponentialAlgorithm, @unpack,
                          AbstractNLSolver, nlsolve_f, issplit,
                          concrete_jac, unwrap_alg, OrdinaryDiffEqCache, _vec, standardtag,
                          isnewton, isnonlinearsolve, _unwrap_val,
                          set_new_W!, set_W_γdt!, alg_difftype, unwrap_cache, diffdir,
                          get_W, isfirstcall, isfirststage, isJcurrent,
                          get_new_W_γdt_cutoff,
                          TryAgain, DIRK, COEFFICIENT_MULTISTEP, NORDSIECK_MULTISTEP, GLM,
                          FastConvergence, Convergence, SlowConvergence,
                          VerySlowConvergence, Divergence, NLStatus, MethodType, constvalue

import OrdinaryDiffEqCore: get_chunksize, resize_J_W!, resize_nlsolver!, alg_autodiff, _get_fwd_tag

using ConstructionBase

import DifferentiationInterface as DI

using FastBroadcast: @..

using ConcreteStructs: @concrete

@static if isdefined(DiffEqBase, :OrdinaryDiffEqTag)
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
