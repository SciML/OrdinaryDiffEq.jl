module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, spzeros, SparseMatrixCSC, AbstractSparseMatrix

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true

# Override the sparse array manipulation functions
OrdinaryDiffEqDifferentiation.nonzeros(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.spzeros(T::Type, m::Integer, n::Integer) = spzeros(T, m, n)

# Helper functions for accessing sparse matrix internals
OrdinaryDiffEqDifferentiation.get_nzval(A::SparseMatrixCSC) = A.nzval
OrdinaryDiffEqDifferentiation.set_all_nzval!(A::SparseMatrixCSC, val) = (A.nzval .= val; A)

end
