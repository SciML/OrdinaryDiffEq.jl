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

# One search of MC21 (CSparse `cs_augment`, iterative): look for an augmenting path from
# the unmatched column `k` through rows `rows[colptr[j]:(limit[j] - 1)]` of each column `j`,
# and flip the matching along it if one is found.
function augment!(matchcol, k, colptr, rows, limit, cheap, visited, search, js, is, ps)
    found = false
    head = 1
    js[1] = k
    i = 0
    while head > 0
        j = js[head]
        if visited[j] != search
            visited[j] = search
            # Cheap assignment: rows before `cheap[j]` are all matched already.
            p = cheap[j]
            while p < limit[j] && !found
                i = rows[p]
                found = matchcol[i] == 0
                p += 1
            end
            cheap[j] = p
            if found
                is[head] = i
                break
            end
            ps[head] = colptr[j]
        end
        p = ps[head]
        while p < limit[j]
            i = rows[p]
            if visited[matchcol[i]] != search
                ps[head] = p + 1
                is[head] = i
                head += 1
                js[head] = matchcol[i]
                break
            end
            p += 1
        end
        p == limit[j] && (head -= 1)
    end
    if found
        for h in head:-1:1
            matchcol[is[h]] = js[h]
        end
    end
    return found
end

function OrdinaryDiffEqDifferentiation.row_matching(
        pattern::SparseMatrixCSC, preferred = nothing
    )
    n = size(pattern, 2)
    size(pattern, 1) == n || return nothing
    colptr, rows, split = candidate_rows(pattern, preferred)
    matchcol = zeros(Int, n)  # column matched to each row, 0 if none
    matched = falses(n)
    cheap = colptr[1:n]
    visited = zeros(Int, n)
    js, is, ps = zeros(Int, n), zeros(Int, n), zeros(Int, n)
    search = 0
    # Match within the preferred entries first, then complete the matching over all of
    # `pattern`. The second pass only reroutes a preferred match on an augmenting path.
    for limit in (split, view(colptr, 2:(n + 1)))
        for k in 1:n
            matched[k] && continue
            search += 1
            matched[k] = augment!(
                matchcol, k, colptr, rows, limit, cheap, visited, search, js, is, ps
            )
        end
    end
    all(matched) || return nothing
    rowof = zeros(Int, n)
    for i in 1:n
        rowof[matchcol[i]] = i
    end
    return rowof
end

end
