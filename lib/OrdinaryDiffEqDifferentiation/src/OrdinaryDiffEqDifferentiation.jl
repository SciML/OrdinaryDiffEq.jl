module OrdinaryDiffEqDifferentiation

import ADTypes
import ADTypes: AutoFiniteDiff, AutoForwardDiff, AbstractADType, AutoSparse

import ForwardDiff, FiniteDiff
import ForwardDiff.Dual
import LinearSolve
import LinearSolve: OperatorAssumptions, LinearVerbosity
import FunctionWrappersWrappers
import DiffEqBase

import LinearAlgebra
import LinearAlgebra: Diagonal, I, UniformScaling, diagind, opnorm
import LinearAlgebra: LowerTriangular
import ArrayInterface

import StaticArraysCore: StaticArray, StaticMatrix

# `_vec`, `_unwrap_val`, `UDerivativeWrapper`, `UJacobianWrapper` are owned by
# SciMLBase (OrdinaryDiffEqCore/DiffEqBase only re-export them), so import them
# from the owner to satisfy `all_explicit_imports_via_owners`.
using SciMLBase: UJacobianWrapper, UDerivativeWrapper, _vec, _unwrap_val
import SciMLBase: SciMLBase, @set, DEIntegrator, ODEFunction, SplitFunction, DAEFunction, islinear, remake, solve!, isconstant
import SciMLOperators: SciMLOperators, update_coefficients, update_coefficients!, MatrixOperator, AbstractSciMLOperator, ScalarOperator
import SparseMatrixColorings: ConstantColoringAlgorithm, GreedyColoringAlgorithm, ColoringProblem,
    ncolors, column_colors, coloring, sparsity_pattern
import OrdinaryDiffEqCore
using OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqAdaptiveImplicitAlgorithm,
    DAEAlgorithm,
    OrdinaryDiffEqImplicitAlgorithm, CompositeAlgorithm,
    OrdinaryDiffEqExponentialAlgorithm,
    OrdinaryDiffEqAdaptiveExponentialAlgorithm,
    AbstractNLSolver, nlsolve_f, issplit,
    concrete_jac, unwrap_alg, OrdinaryDiffEqCache,
    isnewton,
    set_new_W!, set_W_γdt!, diffdir,
    get_W, isfirstcall, isfirststage, isJcurrent,
    get_new_W_γdt_cutoff, isWmethod,
    TryAgain,
    Divergence, constvalue, @SciMLMessage

import OrdinaryDiffEqCore: get_chunksize, resize_J_W!, alg_autodiff

import ConstructionBase

import DifferentiationInterface as DI

using FastBroadcast: @..

using ConcreteStructs: @concrete

import DiffEqBase: OrdinaryDiffEqTag

# Functions for sparse array handling - will be overloaded by extension
# Default implementations return false/error for non-sparse types
is_sparse(::Any) = false
is_sparse_csc(::Any) = false

# These will error if called without the extension, but should never be called
# on non-sparse types due to the is_sparse checks
function nonzeros end
function spzeros end
function get_nzval end
function set_all_nzval! end

# Provide error messages if these are called without extension
nonzeros(A) = error("SparseArrays extension not loaded. Please load SparseArrays to use sparse matrix functionality.")
spzeros(args...) = error("SparseArrays extension not loaded. Please load SparseArrays to use sparse matrix functionality.")
get_nzval(A) = error("SparseArrays extension not loaded. Please load SparseArrays to use sparse matrix functionality.")
set_all_nzval!(A, val) = error("SparseArrays extension not loaded. Please load SparseArrays to use sparse matrix functionality.")

include("alg_utils.jl")
include("linsolve_utils.jl")
include("derivative_utils.jl")
include("derivative_wrappers.jl")
include("operators.jl")

# Documented solver-author extension API: the Jacobian/W-matrix and
# differentiation-config hooks that downstream OrdinaryDiffEq solver sublibraries
# build on. Codegen/perf internals are intentionally left non-public.
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Expr(
            :public,
            :build_J_W, :build_grad_config, :build_jac_config, :build_uf,
            :calc_J, :calc_J!, :calc_tderivative, :calc_tderivative!,
            :calc_rosenbrock_differentiation, :calc_rosenbrock_differentiation!,
            :jacobian!, :jacobian2W!, :update_W!,
            :resize_grad_config!, :resize_jac_config!,
            :dolinsolve, :wrapprecs, :is_always_new, :islinearfunction, :issuccess_W
        )
    )
end

end
