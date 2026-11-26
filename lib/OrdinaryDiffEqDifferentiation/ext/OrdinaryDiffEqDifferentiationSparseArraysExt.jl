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

"""
    structural_matching(A::SparseMatrixCSC) -> Vector{Int}

A maximum matching of the bipartite column-row graph of `A`'s stored pattern, found by
Hopcroft–Karp (SIAM J. Comput. 2(4), 1973) after a greedy start that matches each column
to its stored diagonal entry where it can. Entry `j` is the index into `nonzeros(A)` of
the entry matched to column `j`, or `0` if column `j` is unmatched.
"""
function structural_matching(A::SparseMatrixCSC)
    m, n = size(A)
    rows = rowvals(A)
    col_edge = zeros(Int, n)
    row_col = zeros(Int, m)
    for pass in 1:2, j in 1:n
        col_edge[j] == 0 || continue
        for k in nzrange(A, j)
            i = rows[k]
            if row_col[i] == 0 && (pass == 2 || i == j)
                col_edge[j] = k
                row_col[i] = j
                break
            end
        end
    end
    layer = zeros(Int, n)
    queue = zeros(Int, n)
    stack = zeros(Int, n)
    next_edge = zeros(Int, n)
    while true
        # Breadth-first layering from the free columns, through matched rows.
        qtail = 0
        for j in 1:n
            if col_edge[j] == 0
                layer[j] = 1
                queue[qtail += 1] = j
            else
                layer[j] = 0
            end
        end
        qtail == 0 && break
        reachable = false
        qhead = 0
        while qhead < qtail
            j = queue[qhead += 1]
            for k in nzrange(A, j)
                j2 = row_col[rows[k]]
                if j2 == 0
                    reachable = true
                elseif layer[j2] == 0
                    layer[j2] = layer[j] + 1
                    queue[qtail += 1] = j2
                end
            end
        end
        reachable || break
        # Depth-first search for vertex-disjoint shortest augmenting paths along the layers.
        for j in 1:n
            next_edge[j] = first(nzrange(A, j))
        end
        for j0 in 1:n
            (col_edge[j0] == 0 && layer[j0] == 1) || continue
            depth = 1
            stack[1] = j0
            while depth > 0
                j = stack[depth]
                if next_edge[j] > last(nzrange(A, j))
                    layer[j] = 0
                    depth -= 1
                    continue
                end
                k = next_edge[j]
                next_edge[j] += 1
                j2 = row_col[rows[k]]
                if j2 == 0
                    for d in depth:-1:1
                        jd = stack[d]
                        kd = next_edge[jd] - 1
                        col_edge[jd] = kd
                        row_col[rows[kd]] = jd
                    end
                    for d in 1:depth
                        layer[stack[d]] = 0
                    end
                    break
                elseif layer[j2] == layer[j] + 1
                    stack[depth += 1] = j2
                end
            end
        end
    end
    return col_edge
end

function OrdinaryDiffEqDifferentiation.dae_W_placeholder(J::SparseMatrixCSC)
    W = copy(J)
    vals = nonzeros(W)
    fill!(vals, zero(eltype(W)))
    for k in structural_matching(W)
        k == 0 || (vals[k] = one(eltype(W)))
    end
    return W
end

end
