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
    # A known pattern is fixed by the caller and was checked against `J` when it was built;
    # `NoSparsityDetector` detects nothing, its all-true pattern is only a fallback.
    detector isa Union{
        OrdinaryDiffEqDifferentiation.ADTypes.KnownJacobianSparsityDetector,
        OrdinaryDiffEqDifferentiation.ADTypes.NoSparsityDetector,
    } && return nothing
    # Detectors may return a dense Boolean matrix.
    # Tracing can fail on a residual that the prototype-based AD handles (for example a
    # comparison on a global tracer); validation is then skipped for this refresh.
    detected = try
        OrdinaryDiffEqDifferentiation.ADTypes.jacobian_sparsity(f, fx, x, detector)
    catch
        return nothing
    end
    detected isa SparseMatrixCSC || (detected = SparseMatrixCSC(detected))
    rows = rowvals(J)
    detected_rows = rowvals(detected)
    detected_vals = nonzeros(detected)
    for j in 1:size(detected, 2), idx in nzrange(detected, j)
        # A sparse Boolean pattern can store `false` entries, which detect nothing.
        iszero(detected_vals[idx]) && continue
        i = detected_rows[idx]
        stored = nzrange(J, j)
        k = searchsortedfirst(rows, i, first(stored), last(stored), Base.Order.Forward)
        k <= last(stored) && rows[k] == i && continue
        throw(
            DimensionMismatch(
                "The detected sparsity pattern of the DAE Jacobian includes entry ($i, $j), \
                which the stored Jacobian pattern, fixed when the integrator was built, \
                does not hold. Either `jac_prototype` does not cover the pattern, or the \
                pattern changed during the solve, which is not supported: \
                pass a `KnownJacobianSparsityDetector` with a `jac_prototype` holding the union of all the patterns, or use a dense `autodiff`."
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
