using OrdinaryDiffEqDifferentiation: row_matching, dae_row_matching
using OrdinaryDiffEqBDF, SparseArrays, Test
using ADTypes: AutoForwardDiff, AutoSparse, KnownJacobianSparsityDetector
using SciMLBase: DAEFunction

pattern(n, entries) = sparse(first.(entries), last.(entries), true, n, n)
# `rowof[j]` is the row matched to column `j`; a matching puts a structural nonzero there.
is_matching(P, rowof) = sort(rowof) == 1:size(P, 1) && all(j -> P[rowof[j], j], 1:size(P, 2))

@testset "row_matching" begin
    # A full diagonal matches to the identity even when other perfect matchings exist.
    P = pattern(3, [(1, 1), (2, 2), (3, 3), (1, 2), (2, 1), (3, 1)])
    @test row_matching(P) == 1:3

    # Rows rotated against columns.
    n = 6
    P = pattern(n, [(mod1(j + 1, n), j) for j in 1:n])
    rowof = row_matching(P)
    @test rowof == [mod1(j + 1, n) for j in 1:n]

    # Reaching a perfect matching needs an augmenting path through an earlier choice.
    P = pattern(3, [(1, 1), (2, 1), (1, 2), (3, 3), (3, 2)])
    @test is_matching(P, row_matching(P))

    # Structurally singular: both columns have only row 1.
    @test row_matching(pattern(2, [(1, 1), (1, 2)])) === nothing

    # Preferred entries are matched first: every diagonal is in the full pattern, but the
    # preferred (∂F/∂du) pattern pairs each column with the next row.
    n = 5
    preferred = pattern(n, [(j + 1, j) for j in 1:(n - 1)])
    P = preferred .| pattern(n, [[(j, j) for j in 1:n]; (1, n)])
    rowof = row_matching(P, preferred)
    @test rowof == [2:n; 1]
    @test is_matching(P, rowof)
end

# Residual rows rotated against the unknowns, so ∂F/∂(du) pairs unknown i with row i + 1,
# while ∂F/∂u alone has a full diagonal.
function rotated_chain!(res, du, u, p, t)
    n = length(u)
    for i in 1:(n - 1)
        res[i + 1] = du[i] + 2u[i] - u[i + 1]
    end
    res[1] = u[n] - u[1] / 2
    return nothing
end

function well_ordered!(res, du, u, p, t)
    res[1] = du[1] + u[1] - u[2]
    res[2] = du[2] + u[2] - u[3]
    res[3] = u[3] - u[1] / 2
    return nothing
end

@testset "dae_row_matching" begin
    n = 8
    u = ones(n)
    union_pattern = pattern(
        n, [[(i + 1, i) for i in 1:(n - 1)]; [(i + 1, i + 1) for i in 1:(n - 1)]; (1, n); (1, 1)]
    )
    known = AutoSparse(
        AutoForwardDiff(); sparsity_detector = KnownJacobianSparsityDetector(union_pattern)
    )
    rotated = [2:n; 1]
    # From the autodiff backend's sparsity detector, and from the problem's
    # `jac_prototype`, without an AutoSparse backend.
    @test dae_row_matching(DAEFunction(rotated_chain!), DFBDF(; autodiff = known), u, nothing, 0.0) == rotated
    f = DAEFunction(rotated_chain!; jac_prototype = float.(union_pattern))
    @test dae_row_matching(f, DFBDF(), u, nothing, 0.0) == rotated

    # No pattern available: no matching.
    @test dae_row_matching(DAEFunction(rotated_chain!), DFBDF(), u, nothing, 0.0) === nothing

    # An identity matching is reported as no matching, so a well-ordered DAE is unchanged.
    f = DAEFunction(well_ordered!; jac_prototype = sparse(float.([1 1 0; 0 1 1; 1 0 1])))
    @test dae_row_matching(f, DFBDF(), ones(3), nothing, 0.0) === nothing
end
