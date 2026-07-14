using OrdinaryDiffEqDifferentiation
using SparseArrays
using LinearAlgebra
using Test

# `is_sparse` is true for every AbstractSparseMatrix, so the nzval helpers guarded by
# it must work for every AbstractSparseMatrix as well -- not just SparseMatrixCSC.
# GPU sparse matrices (CuSparseMatrixCSC/CSR/...) are AbstractSparseMatrix with
# non-CSC storage and used to fall through to the "SparseArrays extension not loaded"
# error, even with SparseArrays loaded. The real cuSPARSE matrices are exercised on
# GPU CI (test/gpu/sparse_default_linsolve_tests.jl); here we stand in for them with a
# minimal non-CSC sparse matrix that only provides the AbstractSparseMatrix interface.
struct COOMatrix{Tv, Ti} <: AbstractSparseMatrix{Tv, Ti}
    m::Int
    n::Int
    rows::Vector{Ti}
    cols::Vector{Ti}
    vals::Vector{Tv}
end
Base.size(A::COOMatrix) = (A.m, A.n)
SparseArrays.nonzeros(A::COOMatrix) = A.vals

A_csc = sparse([1.0 0.0; 0.5 2.0])
A_coo = COOMatrix(2, 2, [1, 2, 2], [1, 1, 2], [1.0, 0.5, 2.0])

@test OrdinaryDiffEqDifferentiation.is_sparse(A_csc)
@test OrdinaryDiffEqDifferentiation.is_sparse(A_coo)

@test OrdinaryDiffEqDifferentiation.get_nzval(A_csc) == [1.0, 0.5, 2.0]
@test OrdinaryDiffEqDifferentiation.get_nzval(A_coo) == [1.0, 0.5, 2.0]

@test OrdinaryDiffEqDifferentiation.set_all_nzval!(A_csc, 7.0) === A_csc
@test all(==(7.0), nonzeros(A_csc))

@test OrdinaryDiffEqDifferentiation.set_all_nzval!(A_coo, 7.0) === A_coo
@test all(==(7.0), nonzeros(A_coo))
