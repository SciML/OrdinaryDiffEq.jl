module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, spzeros, SparseMatrixCSC, AbstractSparseMatrix, nnz, sparse

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true

# `nnz` counts stored entries, which is exactly what the colouring consumes, so a prototype
# whose values are all zero but whose structure is real (a `Tridiagonal` of zeros, a
# `SparseMatrixCSC` with stored zeros) is correctly accepted.
function OrdinaryDiffEqDifferentiation._declares_no_nonzeros(S::AbstractSparseMatrix)
    return nnz(S) == 0
end
function OrdinaryDiffEqDifferentiation._declares_no_nonzeros(S::AbstractMatrix)
    return nnz(sparse(S)) == 0
end

# Override the sparse array manipulation functions
OrdinaryDiffEqDifferentiation.nonzeros(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.spzeros(T::Type, m::Integer, n::Integer) = spzeros(T, m, n)

# Helper functions for accessing sparse matrix internals
OrdinaryDiffEqDifferentiation.get_nzval(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.set_all_nzval!(A::AbstractSparseMatrix, val) = (nonzeros(A) .= val; A)

end
