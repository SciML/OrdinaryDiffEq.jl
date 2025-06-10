using Test
using OrdinaryDiffEqRKIP: RKIP
using SciMLBase: SplitODEProblem, SplitFunction, solve
using SciMLOperators: ScalarOperator, DiagonalOperator, MatrixOperator, AbstractSciMLOperator

using StaticArrays

import LinearAlgebra

LinearAlgebra.exp(d::ScalarOperator, t) = ScalarOperator(exp(t * d.val)) # Temporary Fix for missing exp dispatch
LinearAlgebra.exp(d::MatrixOperator, t) = MatrixOperator(exp(t * d.A))

function test(A_prototype, u_prototype, iip; use_ldiv=false)
    for reltol in [1e-3, 1e-6, 1e-8] # fail for 1e-9, probably due to floating point error
        μ = 1.05
        analytic = (u0, _, t) -> u0 .* exp(2μ * t)

        for μₐ in [2μ, μ, 0.0]

            Â = A_prototype(μₐ)
            f = (du, u, _, _) -> du .= (2μ - μₐ) .* u

            if !iip
                f = (u, _, _) -> (2μ - μₐ) .* u
            end

            u0 = u_prototype(1.0)
            t = (0.0, 1.0)

            splfc = SplitFunction{iip}(Â, f; analytic=analytic)
            spltode = SplitODEProblem(splfc, u0, t)
            sol = solve(spltode, RKIP(; use_ldiv=use_ldiv); reltol=reltol)
            @test isapprox(sol(t[end]), splfc.analytic(u0, nothing, t[end]); rtol=reltol)
        end
    end
end

@testset "In-Place ScalarOperator Vector Adaptative Ldiv" begin
    test(μ -> ScalarOperator(μ), u0 -> [u0], true; use_ldiv=true)
end

@testset "In-Place DiagonalOperator Vector Adaptative Ldiv" begin
    test(μ -> DiagonalOperator([μ]), u0 -> [u0], true; use_ldiv=true)
end


@testset "In-Place ScalarOperator Vector Adaptative" begin
    test(μ -> ScalarOperator(μ), u0 -> [u0], true)
end

@testset "In-Place DiagonalOperator Vector Adaptative" begin
    test(μ -> DiagonalOperator([μ]), u0 -> [u0], true)
end

@testset "In-Place MatrixOperator Vector Float Adaptative" begin
    test(μ -> MatrixOperator([μ;;]), u0 -> [u0], true)
end

@testset "Out-of-place Scalar Operator Float Adaptative" begin
    test(μ -> ScalarOperator(μ), u0 -> u0, false)
end

@testset "Out-of-place Scalar Operator SVector Adaptative" begin
    test(μ -> ScalarOperator(μ), u0 -> SVector(u0), false)
end

# Failing Test due to missing dispatch in SciMLOperator


# fails as calling a MatrixOperator on  a StaticVector change its type to Vector, causing a type instability

# @testset "Out-of-place Diagonal 1x1 Operator SVector Adaptative" begin
#     test((μ) -> DiagonalOperator([μ]), u0 -> SVector(u0), false)
# end

# @testset "Out-of-place Matrix 1x1 Operator SVector Adaptative" begin
#     test((μ) -> MatrixOperator([μ;;]), u0 -> SVector(u0), false)
# end

# # Fails for a strange reason as calling ldiv!(a, B, c) on a MatrixOperator dispatch toward an inexisting method 

# @testset "In-Place MatrixOperator Vector Float Adaptative Ldiv" begin
#     test(μ -> MatrixOperator([μ;;]), u0 -> [u0], true; use_ldiv=true)
# end 
