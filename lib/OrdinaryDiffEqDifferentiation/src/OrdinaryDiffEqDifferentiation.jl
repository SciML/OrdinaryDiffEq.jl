module OrdinaryDiffEqDifferentiation

import ADTypes: AutoFiniteDiff, AutoForwardDiff

import SparseDiffTools: SparseDiffTools, matrix_colors, forwarddiff_color_jacobian!,
                forwarddiff_color_jacobian, ForwardColorJacCache,
                default_chunk_size, getsize, JacVec

import ForwardDiff, FiniteDiff
import ForwardDiff.Dual
import LinearSolve
import LinearSolve: OperatorAssumptions
import FunctionWrappersWrappers
using DiffEqBase

import LinearAlgebra
import LinearAlgebra: Diagonal, I, UniformScaling, diagind, mul!, axpby!, opnorm
import SparseArrays: SparseMatrixCSC, AbstractSparseMatrix, nonzeros

import InteractiveUtils
import ArrayInterface

import StaticArrayInterface
import StaticArrays: SArray, MVector, SVector, @SVector, StaticArray, MMatrix, SA,
                     StaticMatrix

using DiffEqBase: TimeGradientWrapper,
                  UJacobianWrapper, TimeDerivativeWrapper,
                  UDerivativeWrapper
using SciMLBase: AbstractSciMLOperator
using OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm, DAEAlgorithm,
OrdinaryDiffEqImplicitAlgorithm, CompositeAlgorithm, OrdinaryDiffEqExponentialAlgorithm,
OrdinaryDiffEqAdaptiveExponentialAlgorithm, @unpack, AbstractNLSolver, nlsolve_f, issplit,
concrete_jac, unwrap_alg, OrdinaryDiffEqCache, _vec, standardtag, isnewton, _unwrap_val,
set_new_W!, set_W_γdt!, alg_difftype, unwrap_cache, diffdir, get_W, isfirstcall, isfirststage, isJcurrent, get_new_W_γdt_cutoff,
TryAgain, DIRK, COEFFICIENT_MULTISTEP, NORDSIECK_MULTISTEP, GLM, FastConvergence, Convergence, SlowConvergence, VerySlowConvergence, Divergence, NLStatus, MethodType

import OrdinaryDiffEq: get_chunksize

import OrdinaryDiffEq: alg_autodiff

using FastBroadcast: @..

@static if isdefined(DiffEqBase, :OrdinaryDiffEqTag)
    import DiffEqBase: OrdinaryDiffEqTag
else
    struct OrdinaryDiffEqTag end
end

include("alg_utils.jl")
include("linsolve_utils.jl")
include("derivative_utils.jl")
include("derivative_wrappers.jl")

end
