module OrdinaryDiffEqCoreSparseArraysExt

using SparseArrays: SparseMatrixCSC
import OrdinaryDiffEqCore: _isdiag, find_algebraic_vars_eqs

# Efficient O(nnz) isdiag check for sparse matrices.
# Standard isdiag is O(n²) which is prohibitively slow for large sparse matrices.
"""
    _isdiag(A::SparseMatrixCSC)

Check if a sparse matrix is diagonal in O(nnz) time by traversing the CSC structure directly.
Returns `true` if all non-zero elements are on the diagonal.
"""
function _isdiag(A::SparseMatrixCSC)
    m, n = size(A)
    m != n && return false
    @inbounds for j in 1:n
        for k in A.colptr[j]:(A.colptr[j + 1] - 1)
            A.rowval[k] != j && return false
        end
    end
    return true
end

"""
    find_algebraic_vars_eqs(M::SparseMatrixCSC)

O(nnz) detection of algebraic variables (zero columns) and equations (zero rows).
"""
function find_algebraic_vars_eqs(M::SparseMatrixCSC)
    n_cols = size(M, 2)
    n_rows = size(M, 1)

    algebraic_vars = fill(true, n_cols)
    algebraic_eqs = fill(true, n_rows)

    @inbounds for j in 1:n_cols
        for idx in M.colptr[j]:(M.colptr[j + 1] - 1)
            if !iszero(M.nzval[idx])
                algebraic_vars[j] = false
                algebraic_eqs[M.rowval[idx]] = false
            end
        end
    end

    return algebraic_vars, algebraic_eqs
end

end
