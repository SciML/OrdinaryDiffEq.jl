module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import ADTypes
import DifferentiationInterface as DI
import OrdinaryDiffEqCore: alg_autodiff
import SparseArrays
import SparseArrays: nonzeros, nzrange, rowvals, spzeros, SparseMatrixCSC, AbstractSparseMatrix

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true
OrdinaryDiffEqDifferentiation.jac_from_pattern(::Type{T}, pattern::SparseMatrixCSC) where {T} =
    SparseMatrixCSC{T}(pattern)

function _missing_entry_error(i, j, reason)
    return DimensionMismatch(
        "$reason entry ($i, $j), which the stored Jacobian pattern, fixed when the \
        integrator was built, does not hold. Either `jac_prototype` does not cover the \
        pattern, or the pattern changed during the solve, which is not supported: pass a \
        `KnownJacobianSparsityDetector` with a `jac_prototype` holding the union of all \
        the patterns, or use a dense `autodiff`."
    )
end

# Whether column `j` of `J` stores row `i`.
function _is_stored(J::SparseMatrixCSC, i, j)
    rows = rowvals(J)
    stored = nzrange(J, j)
    k = searchsortedfirst(rows, i, first(stored), last(stored), Base.Order.Forward)
    return k <= last(stored) && rows[k] == i
end

# The detector cannot trace `f` here, so differentiate it one independent column at a time
# with reused buffers. Only nonzeros outside `J` are errors, since entries that are zero at
# this point add nothing; the values for `J`'s stored entries are staged and copied in once
# every column has passed, so a failure leaves `J` untouched.
function _fill_from_columns!(ad, f, fx, x, J)
    dense = ADTypes.dense_ad(ad)
    tangent = zero(x)
    column = similar(fx)
    prep = try
        DI.prepare_pushforward(f, fx, dense, x, (tangent,); strict = Val(false))
    catch err
        err isa InterruptException && rethrow()
        throw(_unvalidatable_error())
    end
    rows = rowvals(J)
    staged = similar(nonzeros(J))
    for j in axes(J, 2)
        tangent[j] = one(eltype(tangent))
        try
            DI.pushforward!(f, fx, (column,), prep, dense, x, (tangent,))
        catch err
            err isa InterruptException && rethrow()
            throw(_unvalidatable_error())
        finally
            tangent[j] = zero(eltype(tangent))
        end
        for i in eachindex(column)
            iszero(column[i]) || _is_stored(J, i, j) ||
                throw(
                _missing_entry_error(
                    i, j,
                    "The detector cannot trace the residual at this point, and uncompressed \
                    differentiation finds a nonzero derivative at"
                )
            )
        end
        for idx in nzrange(J, j)
            staged[idx] = column[rows[idx]]
        end
    end
    copyto!(nonzeros(J), staged)
    return true
end

_unvalidatable_error() = DimensionMismatch(
    "The automatically detected sparsity pattern of the DAE Jacobian cannot be validated \
    here because the detector cannot trace the residual at this point. Supply a \
    `KnownJacobianSparsityDetector` holding the union of the patterns, or a `jac_prototype`."
)

function OrdinaryDiffEqDifferentiation.check_pattern_stored(
        alg, f, fx, x, J::SparseMatrixCSC, explicit::Bool
    )
    ad = alg_autodiff(alg)
    ad isa ADTypes.AutoSparse || return false
    detector = ADTypes.sparsity_detector(ad)
    # A known pattern is fixed by the caller and was checked against `J` when it was built;
    # `NoSparsityDetector` detects nothing, its all-true pattern is only a fallback.
    detector isa Union{ADTypes.KnownJacobianSparsityDetector, ADTypes.NoSparsityDetector} &&
        return false
    # Storage holding every entry cannot miss one.
    length(nonzeros(J)) == length(J) && return false
    detected = try
        ADTypes.jacobian_sparsity(f, fx, x, detector)
    catch err
        err isa InterruptException && rethrow()
        # Storage the user supplied is theirs to vouch for. Automatically derived storage
        # must not be used unvalidated: a coloring could fold a missing entry into the
        # wrong one.
        explicit && return false
        return _fill_from_columns!(ad, f, fx, x, J)
    end
    # Detectors may return a dense Boolean matrix.
    detected isa SparseMatrixCSC || (detected = SparseMatrixCSC(detected))
    detected_rows = rowvals(detected)
    detected_vals = nonzeros(detected)
    for j in 1:size(detected, 2), idx in nzrange(detected, j)
        # A sparse Boolean pattern can store `false` entries, which detect nothing.
        iszero(detected_vals[idx]) && continue
        i = detected_rows[idx]
        _is_stored(J, i, j) ||
            throw(_missing_entry_error(i, j, "The detected sparsity pattern of the DAE Jacobian includes"))
    end
    return false
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
