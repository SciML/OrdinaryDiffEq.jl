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

# Position of row `i` in the stored entries of column `j` of `J`, or 0 if it is not stored.
function _stored_index(J::SparseMatrixCSC, i, j)
    rows = rowvals(J)
    stored = nzrange(J, j)
    k = searchsortedfirst(rows, i, first(stored), last(stored), Base.Order.Forward)
    return k <= last(stored) && rows[k] == i ? k : 0
end
_is_stored(J::SparseMatrixCSC, i, j) = _stored_index(J, i, j) != 0

# The detector cannot trace `f` here, so differentiate it one independent column at a time
# with reused buffers. Only nonzeros outside `J` are errors, since entries that are zero at
# this point add nothing; the values for `J`'s stored entries are staged and copied in once
# every column has passed, so a failure leaves `J` untouched.
# `f(out, x)` as a function of the single coordinate `x[j]`, so a finite-difference step is
# scaled by that coordinate alone rather than by the norm of the whole state.
mutable struct _CoordinateFunction{F, X}
    f::F
    x::X
    j::Int
end
function (c::_CoordinateFunction)(out, s)
    original = c.x[c.j]
    c.x[c.j] = s[1]
    try
        c.f(out, c.x)
    finally
        c.x[c.j] = original
    end
    return nothing
end

# Returns `column!(column, j)` writing `∂f/∂x[j]`. Finite differences differentiate the
# one-coordinate function with the backend's own stencil and step; AD backends push the
# basis vector `e_j` forward.
function _column_evaluator(dense::ADTypes.AutoFiniteDiff, f, fx, x)
    # A complex-step stencil perturbs the coordinate by an imaginary amount, so the working
    # state must be able to hold it.
    state = dense.fdjtype isa Val{:complex} ? complex.(x) : copy(x)
    coordinate = _CoordinateFunction(f, state, 1)
    s = [first(x)]
    jac = similar(fx, length(fx), 1)
    prep = DI.prepare_jacobian(coordinate, fx, dense, s; strict = Val(false))
    return function (column, j)
        coordinate.j = j
        s[1] = x[j]
        DI.jacobian!(coordinate, fx, jac, prep, dense, s)
        copyto!(column, jac)
        return column
    end
end
function _column_evaluator(dense, f, fx, x)
    tangent = zero(x)
    prep = DI.prepare_pushforward(f, fx, dense, x, (tangent,); strict = Val(false))
    return function (column, j)
        tangent[j] = one(eltype(tangent))
        try
            DI.pushforward!(f, fx, (column,), prep, dense, x, (tangent,))
        finally
            tangent[j] = zero(eltype(tangent))
        end
        return column
    end
end

# Reverse-mode backends differentiate one output row per sweep, so a column at a time would
# sweep every output once per input. Stream output basis vectors instead, through reused
# buffers, staging the stored values of each row.
function _fill_from_rows!(dense, f, fx, x, J)
    cotangent = zero(fx)
    row = similar(x)
    prep = try
        DI.prepare_pullback(f, fx, dense, x, (cotangent,); strict = Val(false))
    catch err
        err isa InterruptException && rethrow()
        throw(_unvalidatable_error())
    end
    staged = zero(nonzeros(J))
    for i in eachindex(fx)
        cotangent[i] = one(eltype(cotangent))
        try
            DI.pullback!(f, fx, (row,), prep, dense, x, (cotangent,))
        catch err
            err isa InterruptException && rethrow()
            throw(_unvalidatable_error())
        finally
            cotangent[i] = zero(eltype(cotangent))
        end
        for j in eachindex(row)
            iszero(row[j]) && continue
            idx = _stored_index(J, i, j)
            idx == 0 && throw(_missing_fallback_entry(i, j))
            staged[idx] = row[j]
        end
    end
    copyto!(nonzeros(J), staged)
    return true
end

_missing_fallback_entry(i, j) = _missing_entry_error(
    i, j,
    "The detector cannot trace the residual at this point, and uncompressed \
    differentiation finds a nonzero derivative at"
)

function _fill_from_columns!(ad, f, fx, x, J)
    dense = ADTypes.dense_ad(ad)
    ADTypes.mode(dense) isa ADTypes.ReverseMode && return _fill_from_rows!(dense, f, fx, x, J)
    column = similar(fx)
    column! = try
        _column_evaluator(dense, f, fx, x)
    catch err
        err isa InterruptException && rethrow()
        throw(_unvalidatable_error())
    end
    rows = rowvals(J)
    staged = similar(nonzeros(J))
    for j in axes(J, 2)
        try
            column!(column, j)
        catch err
            err isa InterruptException && rethrow()
            throw(_unvalidatable_error())
        end
        for i in eachindex(column)
            iszero(column[i]) || _is_stored(J, i, j) ||
                throw(_missing_fallback_entry(i, j))
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
