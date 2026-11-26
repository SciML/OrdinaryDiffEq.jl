module OrdinaryDiffEqDifferentiationSparseArraysExt

using OrdinaryDiffEqDifferentiation
import SparseArrays
import SparseArrays: nonzeros, nzrange, rowvals, spzeros, SparseMatrixCSC, AbstractSparseMatrix

# Override the sparse checking functions
OrdinaryDiffEqDifferentiation.is_sparse(::AbstractSparseMatrix) = true
OrdinaryDiffEqDifferentiation.is_sparse_csc(::SparseMatrixCSC) = true

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

# Candidate rows of each column in the order they are tried: the nonzero entries of
# `preferred` first, then the rest of `pattern`, each group led by the diagonal.
# `split[j]` ends column `j`'s preferred group.
function candidate_rows(pattern::SparseMatrixCSC, preferred)
    n = size(pattern, 2)
    colptr = Vector{Int}(undef, n + 1)
    rows = Vector{Int}(undef, SparseArrays.nnz(pattern))
    split = Vector{Int}(undef, n)
    ispreferred = falses(size(pattern, 1))
    colptr[1] = 1
    for j in 1:n
        if preferred !== nothing
            for p in nzrange(preferred, j)
                iszero(nonzeros(preferred)[p]) || (ispreferred[rowvals(preferred)[p]] = true)
            end
        end
        column = rowvals(pattern)[nzrange(pattern, j)]
        sort!(column; by = i -> (!ispreferred[i], i != j, i))
        copyto!(rows, colptr[j], column)
        colptr[j + 1] = colptr[j] + length(column)
        split[j] = colptr[j] + count(i -> ispreferred[i], column)
        preferred === nothing || (ispreferred[rowvals(preferred)[nzrange(preferred, j)]] .= false)
    end
    return colptr, rows, split
end

# Kuhn's augmenting-path method with Duff's free-row lookahead (I. S. Duff, "On algorithms
# for obtaining a maximum transversal", ACM TOMS 7(3), 1981), as an explicit-stack
# depth-first search. Column `j` searches its candidates `rows[colptr[j]:(stop[j] - 1)]`.
# `colmatch[j]` is the row matched to column `j` and `rowmatch[i]` the column matched to
# row `i`, 0 when unmatched. A row entered in this search is marked with `stamp`.
function augment_from!(rowmatch, colmatch, root, colptr, rows, stop, mark, stamp, stack, cursor)
    depth = 1
    stack[1] = root
    cursor[1] = 0
    while depth > 0
        j = stack[depth]
        if cursor[depth] == 0
            cursor[depth] = colptr[j]
            for p in colptr[j]:(stop[j] - 1)
                rowmatch[rows[p]] == 0 || continue
                # The path's columns move one row along: each takes the row through which
                # it was left, and the deepest takes the free row.
                i = rows[p]
                for d in depth:-1:1
                    c = stack[d]
                    rowmatch[i] = c
                    i, colmatch[c] = colmatch[c], i
                end
                return true
            end
        end
        p = cursor[depth]
        while p < stop[j] && mark[rows[p]] == stamp
            p += 1
        end
        if p < stop[j]
            i = rows[p]
            mark[i] = stamp
            cursor[depth] = p + 1
            depth += 1
            stack[depth] = rowmatch[i]
            cursor[depth] = 0
        else
            depth -= 1
        end
    end
    return false
end

function OrdinaryDiffEqDifferentiation.row_matching(
        pattern::SparseMatrixCSC, preferred = nothing
    )
    n = size(pattern, 2)
    size(pattern, 1) == n || return nothing
    colptr, rows, split = candidate_rows(pattern, preferred)
    rowmatch, colmatch = zeros(Int, n), zeros(Int, n)
    mark, stack, cursor = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    stamp = 0
    # Match within the preferred entries first, then complete the matching over all of
    # `pattern`. An augmenting path keeps every matched column matched, so the second pass
    # only reroutes a preferred match when that is needed to match another column.
    for stop in (split, view(colptr, 2:(n + 1)))
        for k in 1:n
            colmatch[k] == 0 || continue
            stamp += 1
            augment_from!(rowmatch, colmatch, k, colptr, rows, stop, mark, stamp, stack, cursor)
        end
    end
    any(iszero, colmatch) && return nothing
    return colmatch
end

end
