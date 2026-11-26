using ADTypes: AutoFiniteDiff
using CUDA: CUDA, CuArray
using CUDA.CUSPARSE: CuSparseMatrixCSC, CuSparseMatrixCSR
using LinearAlgebra: Diagonal
using OrdinaryDiffEqDefault: ODEFunction, ODEProblem, DefaultImplicitODEAlgorithm,
    solve, successful_retcode
using SparseArrays: sparse
using Test: @test, @testset

CUDA.allowscalar(false)
n = 501
rates = [isodd(i) ? 1.0 : 1000.0 for i in 1:n]
decay!(du, u, p, t) = (@. du = -p * u)
pattern = sparse(Diagonal(ones(n)))
csr = CuSparseMatrixCSR(pattern)
csc = CuSparseMatrixCSC(pattern)

function test_prototype(prototype, rates, slot)
    f = ODEFunction(decay!; jac_prototype = prototype)
    prob = ODEProblem(f, CuArray(ones(length(rates))), (0.0, 0.002), CuArray(rates))
    sol = solve(
        prob, DefaultImplicitODEAlgorithm(autodiff = AutoFiniteDiff());
        abstol = 1.0e-9, reltol = 1.0e-7, maxiters = 1000
    )
    @test successful_retcode(sol)
    @test all(==(slot), sol.alg_choice)
    @test (sol.stats.njacs > 0) == (slot == 5)
    return @test maximum(abs, Array(sol.u[end]) .- exp.(-rates .* 0.002)) < 1.0e-5
end

@testset "cuSPARSE without CUDSS" begin
    test_prototype(csr, rates, 6)
end

using CUDSS: CUDSS

@testset "cuSPARSE with CUDSS" begin
    test_prototype(csr, rates, 5)
    test_prototype(csc, rates, 6)
end
