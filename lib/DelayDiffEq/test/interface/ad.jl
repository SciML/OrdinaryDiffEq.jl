using DelayDiffEq
using OrdinaryDiffEqTsit5

import FiniteDiff
import ForwardDiff

using Test

# Hutchinson's equation
f(u, h, p, t) = p[2] * u * (1 - h(p, t - p[1]) / p[3])
h(p, t) = p[4]

@testset "Gradient" begin
    function test(p)
        prob = DDEProblem(
            f, p[5], h, eltype(p).((0.0, 10.0)), copy(p);
            constant_lags = (p[1],)
        )
        sol = solve(
            prob, MethodOfSteps(Tsit5()); abstol = 1.0e-14, reltol = 1.0e-14,
            saveat = 1.0
        )
        sum(sol)
    end

    # without delay length estimation
    p = [1.5, 1.0, 0.5]
    findiff = FiniteDiff.finite_difference_gradient(p -> test(vcat(1, p, p[end])), p)
    fordiff = ForwardDiff.gradient(p -> test(vcat(1, p, p[end])), p)
    @test maximum(abs.(findiff .- fordiff)) < 1.0e-6

    # with delay length estimation and without discontinuity
    p = [1.0, 1.5, 1.0, 0.5]
    findiff2 = FiniteDiff.finite_difference_gradient(p -> test(vcat(p, p[end])), p)
    fordiff2 = ForwardDiff.gradient(p -> test(vcat(p, p[end])), p)
    @test maximum(abs.(findiff2 .- fordiff2)) < 2

    # with discontinuity and without delay length estimation
    p = [1.5, 1.0, 0.5, 1.0]
    findiff3 = FiniteDiff.finite_difference_gradient(p -> test(vcat(1, p)), p)
    fordiff3 = ForwardDiff.gradient(p -> test(vcat(1, p)), p)
    @test maximum(abs.(findiff3 .- fordiff3)) < 9

    # consistency checks
    @test findiff2[2:end] ≈ findiff
    @test fordiff2[2:end] ≈ fordiff
    @test_broken findiff3[1:(end - 1)] ≈ findiff
    @test_broken fordiff3[1:(end - 1)] ≈ fordiff
end

@testset "Jacobian" begin
    function test(p)
        prob = DDEProblem(
            f, p[5], h, eltype(p).((0.0, 10.0)), copy(p);
            constant_lags = (p[1],)
        )
        sol = solve(
            prob, MethodOfSteps(Tsit5()); abstol = 1.0e-14, reltol = 1.0e-14,
            saveat = 1.0
        )
        sol.u
    end

    # without delay length estimation
    p = [1.5, 1.0, 0.5]
    findiff = FiniteDiff.finite_difference_jacobian(p -> test(vcat(1, p, p[end])), p)
    fordiff = ForwardDiff.jacobian(p -> test(vcat(1, p, p[end])), p)
    @test maximum(abs.(findiff .- fordiff)) < 2.0e-6

    # with delay length estimation and without discontinuity
    p = [1.0, 1.5, 1.0, 0.5]
    findiff2 = FiniteDiff.finite_difference_jacobian(p -> test(vcat(p, p[end])), p)
    fordiff2 = ForwardDiff.jacobian(p -> test(vcat(p, p[end])), p)
    @test maximum(abs.(findiff2 .- fordiff2)) < 3

    # with discontinuity and without delay length estimation
    p = [1.5, 1.0, 0.5, 1.0]
    findiff3 = FiniteDiff.finite_difference_jacobian(p -> test(vcat(1, p)), p)
    fordiff3 = ForwardDiff.jacobian(p -> test(vcat(1, p)), p)
    @test maximum(abs.(findiff3 .- fordiff3)) < 1

    # consistency checks
    @test findiff2[:, 2:end] ≈ findiff
    @test fordiff2[:, 2:end] ≈ fordiff
    @test_broken findiff3[:, 1:(end - 1)] ≈ findiff
    @test_broken fordiff3[:, 1:(end - 1)] ≈ fordiff
end

@testset "Hessian" begin
    function test(p)
        prob = DDEProblem(
            f, p[5], h, eltype(p).((0.0, 10.0)), copy(p);
            constant_lags = (p[1],)
        )
        sol = solve(
            prob, MethodOfSteps(Tsit5()); abstol = 1.0e-14, reltol = 1.0e-14,
            saveat = 1.0
        )
        sum(sol)
    end

    # without delay length estimation and without discontinuity
    p = [1.5, 1.0, 0.5]
    findiff = FiniteDiff.finite_difference_hessian(p -> test(vcat(1, p, p[end])), p)
    fordiff = ForwardDiff.hessian(p -> test(vcat(1, p, p[end])), p)
    @test maximum(abs.(findiff .- fordiff)) < 1.0e-5

    # with delay length estimation and without discontinuity
    p = [1.0, 1.5, 1.0, 0.5]
    findiff2 = FiniteDiff.finite_difference_hessian(p -> test(vcat(p, p[end])), p)
    fordiff2 = ForwardDiff.hessian(p -> test(vcat(p, p[end])), p)
    @test maximum(abs.(findiff2 .- fordiff2)) < 25

    # with discontinuity and without delay length estimation
    p = [1.5, 1.0, 0.5, 1.0]
    findiff3 = FiniteDiff.finite_difference_hessian(p -> test(vcat(1, p)), p)
    fordiff3 = ForwardDiff.hessian(p -> test(vcat(1, p)), p)
    @test maximum(abs.(findiff3 .- fordiff3)) < 1

    # consistency checks
    @test findiff2[2:end, 2:end] ≈ findiff
    @test fordiff2[2:end, 2:end] ≈ fordiff
    @test_broken findiff3[1:(end - 1), 1:(end - 1)] ≈ findiff
    @test_broken fordiff3[1:(end - 1), 1:(end - 1)] ≈ fordiff
end
