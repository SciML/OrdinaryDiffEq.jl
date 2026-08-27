module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, spzeros, SparseMatrixCSC, AbstractSparseMatrix, nnz

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true

function OrdinaryDiffEqDifferentiation._declares_no_nonzeros(S::AbstractSparseMatrix)
    return nnz(S) == 0
end

# Override the sparse array manipulation functions
OrdinaryDiffEqDifferentiation.nonzeros(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.spzeros(T::Type, m::Integer, n::Integer) = spzeros(T, m, n)

# Helper functions for accessing sparse matrix internals
OrdinaryDiffEqDifferentiation.get_nzval(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.set_all_nzval!(A::AbstractSparseMatrix, val) = (nonzeros(A) .= val; A)

end
