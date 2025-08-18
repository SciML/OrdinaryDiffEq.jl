module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, spzeros

# Set the type aliases when extension loads
function __init__()
    OrdinaryDiffEqDifferentiation.SparseMatrixCSC = SparseArrays.SparseMatrixCSC
    OrdinaryDiffEqDifferentiation.AbstractSparseMatrix = SparseArrays.AbstractSparseMatrix
end

# Define functions that were previously imported directly
OrdinaryDiffEqDifferentiation.nonzeros(A::SparseArrays.AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.spzeros(T::Type, m::Integer, n::Integer) = spzeros(T, m, n)

# Helper function to check if a type is sparse
OrdinaryDiffEqDifferentiation.is_sparse_type(::Type{<:SparseArrays.SparseMatrixCSC}) = true
OrdinaryDiffEqDifferentiation.is_sparse_type(::SparseArrays.SparseMatrixCSC) = true
OrdinaryDiffEqDifferentiation.is_sparse_type(::Type{<:SparseArrays.AbstractSparseMatrix}) = true
OrdinaryDiffEqDifferentiation.is_sparse_type(::SparseArrays.AbstractSparseMatrix) = true

# Helper functions for accessing sparse matrix internals
OrdinaryDiffEqDifferentiation.get_nzval(A::SparseArrays.SparseMatrixCSC) = A.nzval
OrdinaryDiffEqDifferentiation.set_all_nzval!(A::SparseArrays.SparseMatrixCSC, val) = (A.nzval .= val; A)

end