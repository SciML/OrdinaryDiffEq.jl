using OrdinaryDiffEqDifferentiation
using SparseArrays
using LinearAlgebra
using SciMLOperators
using SciMLBase
using ADTypes
using Test

f!(du, u, p, t) = (du .= u; nothing)
ad_alg = AutoForwardDiff()

# `findall`/`getindex` on the mass matrix need a concrete matrix; a SciMLOperator
# (e.g. `MatrixOperator`) supports neither directly, so it must be unwrapped first.
@testset "MatrixOperator mass matrix" begin
    M = sparse([1.0 2.0; 0.0 1.0])
    mass_matrix = MatrixOperator(M; update_func = (A, u, p, t) -> M)
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; mass_matrix, jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == Matrix(M)
    @test Matrix(prob.f.jac_prototype) == Matrix(M)
end

@testset "UniformScaling mass matrix (unchanged path)" begin
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == Matrix(jac_prototype)
    @test Matrix(prob.f.jac_prototype) == Matrix(jac_prototype)
end

@testset "Plain sparse mass matrix (unchanged path)" begin
    M = sparse([1.0 2.0; 0.0 1.0])
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; mass_matrix = M, jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == Matrix(M)
    @test Matrix(prob.f.jac_prototype) == Matrix(M)
end

@testset "DiagonalOperator mass matrix" begin
    mass_matrix = DiagonalOperator([1.0, 2.0])
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    sparsity = copy(jac_prototype)
    odef = ODEFunction(f!; mass_matrix, jac_prototype, sparsity)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.sparsity) == [1.0 0.0; 0.0 2.0]
    @test Matrix(prob.f.jac_prototype) == [1.0 0.0; 0.0 2.0]
end

# `sparsity` defaults to (and aliases) `jac_prototype`, which is how #2929's report
# constructs the problem.
@testset "MatrixOperator mass matrix, sparsity defaulted" begin
    M = sparse([1.0 2.0; 0.0 1.0])
    mass_matrix = MatrixOperator(M; update_func = (A, u, p, t) -> M)
    jac_prototype = sparse([1.0 0.0; 0.0 1.0])
    odef = ODEFunction(f!; mass_matrix, jac_prototype)
    prob = ODEProblem(odef, ones(2), (0.0, 1.0))

    OrdinaryDiffEqDifferentiation.prepare_user_sparsity(ad_alg, prob)
    @test Matrix(prob.f.jac_prototype) == Matrix(M)
end

@testset "a prototype declaring no structural nonzeros is rejected" begin
    # AD writes only the stored structure, so a pattern with nothing stored leaves `J`
    # identically zero and the solve is silently wrong (#4140). What matters is the stored
    # count, not the values: a `Tridiagonal` of zeros and a `SparseMatrixCSC` with stored
    # zeros both declare real structure and must still be accepted.
    ad_alg = AutoForwardDiff()
    prob_with(P) = ODEProblem(ODEFunction(f!; jac_prototype = P), ones(2), (0.0, 1.0))

    for P in (
            zeros(2, 2), spzeros(2, 2), Diagonal(zeros(2)),
            Symmetric(zeros(2, 2)), falses(2, 2),
        )
        @test OrdinaryDiffEqDifferentiation._declares_no_nonzeros(P)
        @test_throws ArgumentError OrdinaryDiffEqDifferentiation.prepare_user_sparsity(
            ad_alg, prob_with(P)
        )
    end

    for P in (
            ones(2, 2), Tridiagonal(zeros(1), zeros(2), zeros(1)),
            SparseMatrixCSC(2, 2, [1, 2, 4], [1, 1, 2], [0.0, 0.0, 0.0]),
        )
        @test !OrdinaryDiffEqDifferentiation._declares_no_nonzeros(P)
        @test OrdinaryDiffEqDifferentiation.prepare_user_sparsity(
            ad_alg, prob_with(P)
        ) !== nothing || true
    end

    # An analytic `jac` supplies the Jacobian, so the pattern is irrelevant there.
    jacf = ODEFunction(f!; jac_prototype = zeros(2, 2), jac = (J, u, p, t) -> (J .= 1))
    @test OrdinaryDiffEqDifferentiation.prepare_user_sparsity(
        ad_alg, ODEProblem(jacf, ones(2), (0.0, 1.0))
    ) !== nothing || true
end

