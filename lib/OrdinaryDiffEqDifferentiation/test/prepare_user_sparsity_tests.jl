using OrdinaryDiffEqDifferentiation
using SparseArrays
using SciMLOperators
using SciMLBase
using ADTypes
using LinearAlgebra
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
