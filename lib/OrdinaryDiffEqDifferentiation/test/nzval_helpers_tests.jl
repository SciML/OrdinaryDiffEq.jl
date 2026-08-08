using OrdinaryDiffEqDifferentiation
using OrdinaryDiffEqRosenbrock
using SparseArrays
using LinearAlgebra
using Test

# `is_sparse` is true for every AbstractSparseMatrix, so the nzval helpers guarded by
# it must work for every AbstractSparseMatrix as well -- not just SparseMatrixCSC.
# GPU sparse matrices (CuSparseMatrixCSC/CSR/...) are AbstractSparseMatrix with
# non-CSC storage and used to fall through to the "SparseArrays extension not loaded"
# error, even with SparseArrays loaded. The real cuSPARSE matrices are exercised on
# GPU CI (test/gpu/sparse_default_linsolve_tests.jl); here we stand in for them with a
# minimal non-CSC sparse matrix that only provides the AbstractSparseMatrix interface.
struct COOMatrix{Tv, Ti} <: AbstractSparseMatrix{Tv, Ti}
    m::Int
    n::Int
    rows::Vector{Ti}
    cols::Vector{Ti}
    vals::Vector{Tv}
end
Base.size(A::COOMatrix) = (A.m, A.n)
SparseArrays.nonzeros(A::COOMatrix) = A.vals

A_csc = sparse([1.0 0.0; 0.5 2.0])
A_coo = COOMatrix(2, 2, [1, 2, 2], [1, 1, 2], [1.0, 0.5, 2.0])

@test OrdinaryDiffEqDifferentiation.is_sparse(A_csc)
@test OrdinaryDiffEqDifferentiation.is_sparse(A_coo)

@test OrdinaryDiffEqDifferentiation.get_nzval(A_csc) == [1.0, 0.5, 2.0]
@test OrdinaryDiffEqDifferentiation.get_nzval(A_coo) == [1.0, 0.5, 2.0]

@test OrdinaryDiffEqDifferentiation.set_all_nzval!(A_csc, 7.0) === A_csc
@test all(==(7.0), nonzeros(A_csc))

@test OrdinaryDiffEqDifferentiation.set_all_nzval!(A_coo, 7.0) === A_coo
@test all(==(7.0), nonzeros(A_coo))

using OrdinaryDiffEqDifferentiation: fill_stored!

@testset "fill_stored!" begin
    A = zeros(2, 2)
    @test fill_stored!(A, 3.0) === A
    @test all(==(3.0), A)

    D = Diagonal(zeros(2))
    @test fill_stored!(D, 3.0) === D
    @test D == Diagonal([3.0, 3.0])

    B = Bidiagonal(zeros(3), zeros(2), :U)
    @test fill_stored!(B, 3.0) === B
    @test B == Bidiagonal([3.0, 3.0, 3.0], [3.0, 3.0], :U)

    T = Tridiagonal(zeros(2), zeros(3), zeros(2))
    @test fill_stored!(T, 3.0) === T
    @test T == Tridiagonal([3.0, 3.0], [3.0, 3.0, 3.0], [3.0, 3.0])

    S = SymTridiagonal(zeros(3), zeros(2))
    @test fill_stored!(S, 3.0) === S
    @test S == SymTridiagonal([3.0, 3.0, 3.0], [3.0, 3.0])
end

# build_J_W used to `fill!` non-sparse jac_prototype caches with one, which
# throws for banded structured matrices whose off-band entries are constrained.
@testset "structured jac_prototype J/W caches initialize with stored ones" begin
    f! = (du, u, p, t) -> du .= u ./ t
    # Zero the off-diagonal stored entries too: they are initialized to one, so
    # a jac! that skips them would leave a wrong Jacobian behind.
    jac! = (J, u, p, t) -> (J[1, 1] = 1 / t; J[2, 2] = 1 / t; J[1, 2] = 0; J[2, 1] = 0; nothing)

    for jac_prototype in (Diagonal(zeros(2)), Tridiagonal(zeros(1), zeros(2), zeros(1)))
        fun = ODEFunction(f!; jac = jac!, jac_prototype = jac_prototype)
        prob = ODEProblem(fun, ones(2), (1.0, 10.0))
        integ = init(prob, Rosenbrock23())
        ones_prototype = fill_stored!(copy(jac_prototype), 1.0)
        @test integ.cache.J isa typeof(jac_prototype)
        @test integ.cache.J == ones_prototype
        @test integ.cache.W == ones_prototype
        sol = solve!(integ)
        @test sol.u[end] ≈ [10.0, 10.0]
    end
end
