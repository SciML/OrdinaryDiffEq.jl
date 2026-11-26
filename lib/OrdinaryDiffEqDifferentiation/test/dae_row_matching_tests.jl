using OrdinaryDiffEqDifferentiation: row_matching, dae_row_matching, supports_row_matching
using OrdinaryDiffEqBDF, SparseArrays, Random, Test
using LinearAlgebra: I
using ADTypes: AutoForwardDiff, AutoSparse, KnownJacobianSparsityDetector
using LinearSolve: KrylovJL_GMRES, KrylovJL_FGMRES, KrylovJL_BICGSTAB, KrylovJL_CG,
    KrylovJL_MINRES
using SciMLBase: DAEFunction

pattern(n, entries) = sparse(first.(entries), last.(entries), true, n, n)
# `rowof[j]` is the row matched to column `j`; a matching puts a structural nonzero there.
is_matching(P, rowof) = sort(rowof) == 1:size(P, 1) && all(j -> P[rowof[j], j], 1:size(P, 2))

# Whether columns `j:n` can be matched to distinct rows not in `used`, by trying every row.
function has_perfect_matching(P, j = 1, used = falses(size(P, 1)))
    j > size(P, 2) && return true
    for i in 1:size(P, 1)
        if P[i, j] && !used[i]
            used[i] = true
            has_perfect_matching(P, j + 1, used) && return true
            used[i] = false
        end
    end
    return false
end

function check_against_exhaustive(P, preferred = nothing)
    rowof = row_matching(sparse(P), preferred === nothing ? nothing : sparse(preferred))
    return has_perfect_matching(P) ? rowof !== nothing && is_matching(P, rowof) :
        rowof === nothing
end

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

    # Matching column 3 reroutes column 1 from its preferred row 1 to row 2, and column 2
    # from its preferred row 2 to row 3.
    P = pattern(3, [(1, 1), (2, 1), (2, 2), (3, 2), (1, 3)])
    preferred = pattern(3, [(1, 1), (2, 2)])
    @test row_matching(P, preferred) == [2, 3, 1]
end

@testset "row_matching against exhaustive search" begin
    # Every pattern up to 3 x 3, with and without its lower triangle preferred.
    for n in 1:3, bits in 0:(2^(n * n) - 1)
        P = reshape([isodd(bits >> k) for k in 0:(n * n - 1)], n, n)
        @test check_against_exhaustive(P)
        @test check_against_exhaustive(P, P .& [i >= j for i in 1:n, j in 1:n])
    end
    rng = Xoshiro(4868)
    for n in 4:7, density in (0.2, 0.35, 0.5), _ in 1:150
        P = rand(rng, n, n) .< density
        @test check_against_exhaustive(P)
        @test check_against_exhaustive(P, P .& (rand(rng, n, n) .< 0.5))
    end
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
    matching = dae_row_matching(DAEFunction(rotated_chain!), DFBDF(; autodiff = known), u, nothing, 0.0)
    @test matching.rowof == rotated
    # Unknown `j < n` is matched to row `j + 1`, whose largest entry is `∂F/∂u = 2`; unknown
    # `n` to row 1, whose largest entry is 1.
    @test matching.rowscale == [fill(0.5, n - 1); 1.0]
    f = DAEFunction(rotated_chain!; jac_prototype = float.(union_pattern))
    @test dae_row_matching(f, DFBDF(), u, nothing, 0.0).rowof == rotated

    # Each row is scaled by its own magnitude, so a row written 1e-14 times smaller gets the
    # same preconditioned residual.
    scaled!(res, du, u, p, t) = (rotated_chain!(res, du, u, p, t); res[1] *= 1.0e-14; nothing)
    scaled = dae_row_matching(DAEFunction(scaled!), DFBDF(; autodiff = known), u, nothing, 0.0)
    @test scaled.rowof == rotated
    @test scaled.rowscale ≈ [fill(0.5, n - 1); 1.0e14]

    # No pattern available, or only a dense one: no matching.
    @test dae_row_matching(DAEFunction(rotated_chain!), DFBDF(), u, nothing, 0.0) === nothing
    dense = AutoSparse(
        AutoForwardDiff(); sparsity_detector = KnownJacobianSparsityDetector(Matrix(union_pattern))
    )
    @test dae_row_matching(DAEFunction(rotated_chain!), DFBDF(; autodiff = dense), u, nothing, 0.0) ===
        nothing

    # An identity matching is reported as no matching, so a well-ordered DAE is unchanged.
    f = DAEFunction(well_ordered!; jac_prototype = sparse(float.([1 1 0; 0 1 1; 1 0 1])))
    @test dae_row_matching(f, DFBDF(), ones(3), nothing, 0.0) === nothing
end

@testset "supports_row_matching" begin
    for linsolve in (KrylovJL_GMRES(), KrylovJL_FGMRES(), KrylovJL_BICGSTAB())
        @test supports_row_matching(linsolve)
    end
    # A row permutation is not a symmetric positive definite preconditioner.
    @test !supports_row_matching(KrylovJL_CG())
    @test !supports_row_matching(KrylovJL_MINRES())
    # A user preconditioner, or one applied with `mul!`, replaces the default one.
    @test !supports_row_matching(KrylovJL_GMRES(; precs = Returns((I, I))))
    @test !supports_row_matching(KrylovJL_GMRES(; ldiv = false))
    @test !supports_row_matching(nothing)
end
