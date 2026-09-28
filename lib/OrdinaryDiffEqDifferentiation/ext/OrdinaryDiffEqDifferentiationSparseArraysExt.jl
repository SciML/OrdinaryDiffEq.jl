module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, nzrange, rowvals, spzeros, SparseMatrixCSC, AbstractSparseMatrix, nnz
import LinearAlgebra: Diagonal, UniformScaling
import SciMLOperators: ScalarOperator

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true

# Override the sparse array manipulation functions
OrdinaryDiffEqDifferentiation.nonzeros(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.spzeros(T::Type, m::Integer, n::Integer) = spzeros(T, m, n)

# Helper functions for accessing sparse matrix internals
OrdinaryDiffEqDifferentiation.get_nzval(A::AbstractSparseMatrix) = nonzeros(A)
OrdinaryDiffEqDifferentiation.set_all_nzval!(A::AbstractSparseMatrix, val) = (nonzeros(A) .= val; A)

# ---------------------------------------------------------------------------
# Pattern-stable W assembly for CPU `SparseMatrixCSC`.
#
# Broadcasting `W = muladd(-M, invdtgamma, J)` into a sparse `W` drops exact
# zeros, so `nnz(W)` fluctuates with the Jacobian's *values* and KLU/UMFPACK
# re-run symbolic analysis nearly every step. These methods keep `W`'s pattern
# fixed instead: the first call widens it to `W ∪ J ∪ nz(M)` (or
# `W ∪ J_u ∪ J_du`), and every call rewrites `nonzeros(W)` over the sorted row
# indices — positionally when `W` and `J` already share a pattern, by a merge
# pass otherwise. Values use the same `muladd` operands as the broadcast path,
# so they are bitwise identical.
# ---------------------------------------------------------------------------

const _ScalarOrDiagonalMassMatrix = Union{UniformScaling, ScalarOperator, Diagonal}

# Every row index stored in `S`'s column `j` must also be stored in `W`'s.
@inline function _column_covers(Wrows, wrange, Srows, srange)
    wp = first(wrange)
    wend = last(wrange)
    @inbounds for sp in srange
        r = Srows[sp]
        while wp <= wend && Wrows[wp] < r
            wp += 1
        end
        (wp <= wend && Wrows[wp] == r) || return false
    end
    return true
end

# Which diagonal slots the mass matrix needs in `W`'s pattern: all of them for
# a scalar `M` (`W = J - λ/dtγ * I` materializes the whole diagonal), only the
# structurally nonzero ones (`nz(M)`) for a `Diagonal`. `parent(::Diagonal)` is
# the backing vector; `diag(::Diagonal)` would allocate a fresh one per call.
_diag_vec_nonzero(d, j) = !iszero(d[j])
_diag_needed(mass_matrix::Diagonal) = Base.Fix1(_diag_vec_nonzero, parent(mass_matrix))
_diag_needed(mass_matrix) = Returns(true)

# Position of `r` in `rows[range]`, or `0` when absent. CSC column indices are
# sorted, so binary search applies.
@inline function _row_pos(rows, r, range)
    lo = first(range)
    hi = last(range)
    @inbounds while lo < hi
        mid = (lo + hi) >>> 1
        if rows[mid] < r
            lo = mid + 1
        else
            hi = mid
        end
    end
    return lo <= hi && rows[lo] == r ? lo : 0
end

@inline _row_in_range(rows, r, range) = _row_pos(rows, r, range) != 0

# Identical nonzero structure: same stored rows and column extents.
function _same_pattern(A::SparseMatrixCSC, B::SparseMatrixCSC)
    size(A) == size(B) && nnz(A) == nnz(B) || return false
    rowvals(A) == rowvals(B) || return false
    @inbounds for j in axes(A, 2)
        length(nzrange(A, j)) == length(nzrange(B, j)) || return false
    end
    return true
end

# Whether `W`'s pattern covers every stored index of each matrix in `srcs` plus
# the diagonal slots `diag_needed` asks for.
function _pattern_covers(W::SparseMatrixCSC, srcs::Tuple, diag_needed)
    Wrows = rowvals(W)
    @inbounds for j in axes(W, 2)
        wrange = nzrange(W, j)
        for S in srcs
            _column_covers(Wrows, wrange, rowvals(S), nzrange(S, j)) ||
                return false
        end
        if j <= size(W, 1) && diag_needed(j) && !_row_in_range(Wrows, j, wrange)
            return false
        end
    end
    return true
end

# Sorted union of the sorted index lists `a` and `b`, appended to `out`.
function _merge_sorted!(out, a, b)
    i = firstindex(a)
    la = lastindex(a)
    j = firstindex(b)
    lb = lastindex(b)
    @inbounds while i <= la || j <= lb
        ra = i <= la ? a[i] : typemax(eltype(a))
        rb = j <= lb ? b[j] : typemax(eltype(b))
        if ra == rb
            push!(out, ra)
            i += 1
            j += 1
        elseif ra < rb
            push!(out, ra)
            i += 1
        else
            push!(out, rb)
            j += 1
        end
    end
    return out
end

# Give `W` the pattern `W ∪ ⋃ srcs ∪ {j : diag_needed(j)}`, preserving its
# stored values. Runs once per solve; afterwards the pattern is fixed.
@noinline function _widen_pattern!(
        W::SparseMatrixCSC{Tv, Ti}, srcs::Tuple, diag_needed
    ) where {Tv, Ti}
    m, n = size(W)
    Wrows = rowvals(W)
    Wnz = nonzeros(W)
    colptr = Vector{Ti}(undef, n + 1)
    colptr[1] = one(Ti)
    rows = Ti[]
    vals = Tv[]
    for j in 1:n
        wrange = nzrange(W, j)
        wcol = @view(Wrows[wrange])
        merged = wcol
        for S in srcs
            merged = _merge_sorted!(
                Ti[], merged, @view(rowvals(S)[nzrange(S, j)])
            )
        end
        if j <= m && diag_needed(j)
            merged = _merge_sorted!(Ti[], merged, Ti[j])
        end
        colptr[j + 1] = colptr[j] + length(merged)
        for r in merged
            push!(rows, r)
            k = searchsortedfirst(wcol, r)
            push!(
                vals,
                k <= length(wcol) && wcol[k] == r ?
                    Wnz[first(wrange) + k - 1] : zero(Tv)
            )
        end
    end
    copyto!(W, SparseMatrixCSC{Tv, Ti}(m, n, colptr, rows, vals))
    return W
end

# Rewrite `nonzeros(W)` as `J - M/dtgamma` over `W`'s fixed pattern by merging
# `J`'s row indices into `W`'s column by column. Slots `J` doesn't store count
# as zero; the diagonal slot adds the mass-matrix contribution.
function _fill_J_minus_M!(
        W::SparseMatrixCSC, mass_matrix::_ScalarOrDiagonalMassMatrix,
        invdtgamma, J::SparseMatrixCSC
    )
    Wrows = rowvals(W)
    Wnz = nonzeros(W)
    Jrows = rowvals(J)
    Jnz = nonzeros(J)
    if mass_matrix isa Diagonal
        d = parent(mass_matrix)
        @inbounds for j in axes(W, 2)
            # Broadcast operand is `-M`: `+0` off the diagonal, `-d[j]` on it.
            aj = -d[j]
            jend = last(nzrange(J, j))
            jp = first(nzrange(J, j))
            for p in nzrange(W, j)
                i = Wrows[p]
                jv = zero(eltype(Jnz))
                if jp <= jend && Jrows[jp] == i
                    jv = Jnz[jp]
                    jp += 1
                end
                Wnz[p] = muladd(i == j ? aj : zero(aj), invdtgamma, jv)
            end
        end
    else
        λ = -OrdinaryDiffEqDifferentiation._scalar_massmatrix_λ(mass_matrix)
        # Scalar `M` copies `J` bitwise, so off-diagonal slots are the bare
        # `jv` (even a stored `-0.0`); only the diagonal takes the shift.
        @inbounds for j in axes(W, 2)
            jend = last(nzrange(J, j))
            jp = first(nzrange(J, j))
            for p in nzrange(W, j)
                i = Wrows[p]
                jv = zero(eltype(Jnz))
                if jp <= jend && Jrows[jp] == i
                    jv = Jnz[jp]
                    jp += 1
                end
                Wnz[p] = i == j ? muladd(λ, invdtgamma, jv) : jv
            end
        end
    end
    return W
end

# `W` and `J` share one pattern: values are a positional rewrite, and each
# column's diagonal slot is met during the walk. Off-diagonal values reproduce
# the caller's path bitwise — a scalar `M` copies `J` verbatim (a stored
# `-0.0` stays), while `muladd(-M, …)` broadcast supplies `+0.0` for a
# `Diagonal`'s off-diagonal slots (`fma(+0.0, c, -0.0) == +0.0`). Both return
# `false` when a needed diagonal slot is missing so the caller can widen the
# pattern first.
function _fill_same_pattern!(
        W::SparseMatrixCSC, mass_matrix, invdtgamma, J::SparseMatrixCSC
    )
    Wrows = rowvals(W)
    Wnz = nonzeros(W)
    Jnz = nonzeros(J)
    m = size(W, 1)
    λ = -OrdinaryDiffEqDifferentiation._scalar_massmatrix_λ(mass_matrix)
    @inbounds for j in axes(W, 2)
        found = false
        for p in nzrange(W, j)
            i = Wrows[p]
            if i == j
                found = true
                Wnz[p] = muladd(λ, invdtgamma, Jnz[p])
            else
                Wnz[p] = Jnz[p]
            end
        end
        j <= m && !found && return false
    end
    return true
end

function _fill_same_pattern!(
        W::SparseMatrixCSC, mass_matrix::Diagonal, invdtgamma, J::SparseMatrixCSC
    )
    Wrows = rowvals(W)
    Wnz = nonzeros(W)
    Jnz = nonzeros(J)
    d = parent(mass_matrix)
    m = size(W, 1)
    @inbounds for j in axes(W, 2)
        aj = j <= m ? -d[j] : zero(eltype(d))
        found = false
        for p in nzrange(W, j)
            i = Wrows[p]
            if i == j
                found = true
                Wnz[p] = muladd(aj, invdtgamma, Jnz[p])
            else
                Wnz[p] = muladd(zero(aj), invdtgamma, Jnz[p])
            end
        end
        j <= m && !found && !iszero(d[j]) && return false
    end
    return true
end

# Rewrite `nonzeros(W)` as `f` of the per-slot values of `S1` and `S2` over
# `W`'s fixed pattern. `W` may alias `S1` (for `A .-= B`): reads precede the
# write at each position, and both walks advance in row order.
function _fill_combine!(
        f, W::SparseMatrixCSC, S1::SparseMatrixCSC, S2::SparseMatrixCSC
    )
    Wrows = rowvals(W)
    Wnz = nonzeros(W)
    r1 = rowvals(S1)
    v1 = nonzeros(S1)
    r2 = rowvals(S2)
    v2 = nonzeros(S2)
    @inbounds for j in axes(W, 2)
        e1 = last(nzrange(S1, j))
        p1 = first(nzrange(S1, j))
        e2 = last(nzrange(S2, j))
        p2 = first(nzrange(S2, j))
        for p in nzrange(W, j)
            i = Wrows[p]
            a = zero(eltype(v1))
            if p1 <= e1 && r1[p1] == i
                a = v1[p1]
                p1 += 1
            end
            b = zero(eltype(v2))
            if p2 <= e2 && r2[p2] == i
                b = v2[p2]
                p2 += 1
            end
            Wnz[p] = f(a, b)
        end
    end
    return W
end

function OrdinaryDiffEqDifferentiation.jacobian2W!(
        W::SparseMatrixCSC, mass_matrix::_ScalarOrDiagonalMassMatrix,
        dtgamma::Number, J::SparseMatrixCSC
    )::Nothing
    iijj = axes(W)
    @boundscheck (iijj == axes(J) && length(iijj) == 2) ||
        OrdinaryDiffEqDifferentiation._throwWJerror(W, J)
    OrdinaryDiffEqDifferentiation._is_scalar_massmatrix(mass_matrix) ||
        @boundscheck axes(mass_matrix) == axes(W) ||
        OrdinaryDiffEqDifferentiation._throwWMerror(W, mass_matrix)
    @inbounds begin
        invdtgamma = inv(dtgamma)
        if _same_pattern(W, J) &&
                _fill_same_pattern!(W, mass_matrix, invdtgamma, J)
            return nothing
        end
        diag_needed = _diag_needed(mass_matrix)
        _pattern_covers(W, (J,), diag_needed) ||
            _widen_pattern!(W, (J,), diag_needed)
        _fill_J_minus_M!(W, mass_matrix, invdtgamma, J)
    end
    return nothing
end

function OrdinaryDiffEqDifferentiation.dae_jacobian2W!(
        W::SparseMatrixCSC, J_u::SparseMatrixCSC,
        J_du::SparseMatrixCSC, cj::Number
    )::Nothing
    @boundscheck axes(W) == axes(J_u) == axes(J_du) ||
        throw(DimensionMismatch("W, J_u, J_du must have matching axes"))
    @inbounds begin
        if _same_pattern(W, J_u) && _same_pattern(W, J_du)
            nonzeros(W) .= muladd.(cj, nonzeros(J_du), nonzeros(J_u))
            return nothing
        end
        _pattern_covers(W, (J_u, J_du), Returns(false)) ||
            _widen_pattern!(W, (J_u, J_du), Returns(false))
        _fill_combine!(W, J_u, J_du) do vu, vdu
            muladd(cj, vdu, vu)
        end
    end
    return nothing
end

# `A .-= B` over a pattern widened once to `A ∪ B`. `A`'s own stored values are
# carried over by `_widen_pattern!`, so the fill can read them back.
function OrdinaryDiffEqDifferentiation._sparse_pattern_stable_sub!(
        A::SparseMatrixCSC, B::SparseMatrixCSC
    )
    @boundscheck axes(A) == axes(B) ||
        throw(DimensionMismatch("A, B must have matching axes"))
    @inbounds begin
        if _same_pattern(A, B)
            nonzeros(A) .-= nonzeros(B)
            return true
        end
        _pattern_covers(A, (B,), Returns(false)) ||
            _widen_pattern!(A, (B,), Returns(false))
        _fill_combine!(-, A, A, B)
    end
    return true
end

end
