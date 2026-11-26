using OrdinaryDiffEqDifferentiation, SparseArrays, LinearAlgebra, Random, Test
using OrdinaryDiffEqDifferentiation: dae_W_placeholder, weight_left, weight_right

const SAExt = Base.get_extension(
    OrdinaryDiffEqDifferentiation, :OrdinaryDiffEqDifferentiationSparseArraysExt
)

function has_perfect_matching(P)
    n = size(P, 2)
    return any(p -> all(j -> P[p[j], j], 1:n), permutations_of(n))
end
permutations_of(n) = n == 1 ? [[1]] :
    [insert!(copy(p), k, n) for p in permutations_of(n - 1) for k in 1:n]

@testset "Structural matching of a sparse pattern" begin
    rng = Xoshiro(4866)
    perfect = 0
    for n in 1:6, density in (0.2, 0.35, 0.6), _ in 1:40
        P = rand(rng, n, n) .< density
        A = sparse(Float64.(P))
        edges = SAExt.structural_matching(A)
        matched = filter(!=(0), edges)
        @test all(j -> edges[j] == 0 || edges[j] in nzrange(A, j), 1:n)
        @test allunique(rowvals(A)[matched])
        W = dae_W_placeholder(A)
        @test W.colptr == A.colptr && rowvals(W) == rowvals(A)
        if has_perfect_matching(P)
            perfect += 1
            @test length(matched) == n
            @test abs(det(Matrix(W))) == 1
        else
            @test length(matched) < n
        end
        if all(i -> P[i, i], 1:n)
            @test Matrix(W) == I
        end
    end
    @test perfect > 50

    # Rerouting: greedy takes row 1 for column 1, and only an augmenting path through
    # column 1 frees it for column 2.
    A = sparse([1, 2, 1], [1, 1, 2], ones(3), 2, 2)
    @test Matrix(dae_W_placeholder(A)) == [0 1; 1 0]
    @test Matrix(dae_W_placeholder(sparse(Float64[1 1; 1 1]))) == I
end

@testset "Residual weights composed with a user preconditioner" begin
    rng = Xoshiro(1)
    n = 5
    w = 1 .+ 9 .* rand(rng, n)
    Pl = lu(rand(rng, n, n) + n * I)
    Pr = lu(UpperTriangular(rand(rng, n, n)) + 2I)
    x = randn(rng, n)
    left = weight_left(Pl, w)
    right = weight_right(Pr, w)
    @test ldiv!(similar(x), left, x) ≈ w .* (Pl \ x)
    @test ldiv!(similar(x), right, x) ≈ Pr \ (x ./ w)
    @test ldiv!(left, copy(x)) ≈ w .* (Pl \ x)
    @test ldiv!(right, copy(x)) ≈ Pr \ (x ./ w)
    @test eltype(left) == eltype(right) == Float64
    @test ldiv!(similar(x), weight_left(I, w), x) ≈ w .* x
    @test weight_right(I, w) == Diagonal(w)
end
