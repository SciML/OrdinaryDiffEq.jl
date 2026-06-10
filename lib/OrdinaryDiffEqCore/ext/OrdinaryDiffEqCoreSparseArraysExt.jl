module OrdinaryDiffEqCoreSparseArraysExt

using SparseArrays: SparseMatrixCSC
import OrdinaryDiffEqCore: _isdiag, _find_large_jac_entries!

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

# only look at nonzero vals
function _find_large_jac_entries!(rows::Set{Int}, cols::Set{Int}, jac::SparseMatrixCSC)
    @inbounds for j in axes(jac, 2)
        for k in jac.colptr[j]:(jac.colptr[j + 1] - 1)
            val = jac.nzval[k]
            if !isfinite(val) || abs(val) > 1e6
                push!(rows, jac.rowval[k])
                push!(cols, j)
            end
        end
    end
end

end
