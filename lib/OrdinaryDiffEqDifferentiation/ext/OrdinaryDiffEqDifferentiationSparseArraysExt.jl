module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, nzrange, rowvals, spzeros, SparseMatrixCSC, AbstractSparseMatrix

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true
OrdinaryDiffEqDifferentiation.jac_from_pattern(::Type{T}, pattern::SparseMatrixCSC) where {T} =
    SparseMatrixCSC{T}(pattern)

function OrdinaryDiffEqDifferentiation.check_pattern_stored(
        alg, f, fx, x, J::SparseMatrixCSC
    )
    ad = OrdinaryDiffEqDifferentiation.alg_autodiff(alg)
    ad isa OrdinaryDiffEqDifferentiation.AutoSparse || return nothing
    detector = OrdinaryDiffEqDifferentiation.ADTypes.sparsity_detector(ad)
    # A known pattern is fixed by the caller and was checked against `J` when it was built.
    detector isa OrdinaryDiffEqDifferentiation.ADTypes.KnownJacobianSparsityDetector &&
        return nothing
    detected = OrdinaryDiffEqDifferentiation.ADTypes.jacobian_sparsity(f, fx, x, detector)
    rows = rowvals(J)
    for j in 1:size(detected, 2), i in rowvals(detected)[nzrange(detected, j)]
        range = nzrange(J, j)
        k = searchsortedfirst(view(rows, range), i)
        k <= length(range) && rows[range[k]] == i && continue
        throw(
            DimensionMismatch(
                "The automatically detected sparsity pattern of the DAE Jacobian grew to \
                include entry ($i, $j), which the stored Jacobian pattern, fixed when the \
                integrator was built, does not hold. The pattern must not change during \
                the solve: pass a `KnownJacobianSparsityDetector` with a `jac_prototype` \
                holding the union of all the patterns, or use a dense `autodiff`."
            )
        )
    end
    return nothing
end

# Override the sparse array manipulation functions
OrdinaryDiffEqDifferentiation.nonzeros(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.spzeros(T::Type, m::Integer, n::Integer) = spzeros(T, m, n)

# Helper functions for accessing sparse matrix internals
OrdinaryDiffEqDifferentiation.get_nzval(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.set_all_nzval!(A::AbstractSparseMatrix, val) = (nonzeros(A) .= val; A)

function OrdinaryDiffEqDifferentiation._update_sparse_diagonal!(
        W::SparseMatrixCSC, λ, invdtgamma, J
    )
    n = size(W, 2)
    size(W, 1) == n || throw(BoundsError(W, (size(W, 1), n)))
    rows = rowvals(W)
    vals = nonzeros(W)
    @inbounds for j in 1:n
        colrange = nzrange(W, j)
        localidx = searchsortedfirst(@view(rows[colrange]), j)
        if localidx <= length(colrange) && rows[first(colrange) + localidx - 1] == j
            pos = first(colrange) + localidx - 1
            vals[pos] = muladd(λ, invdtgamma, vals[pos])
        else
            for i in j:n
                W[i, i] = muladd(λ, invdtgamma, J[i, i])
            end
            return true
        end
    end
    return true
end

end
