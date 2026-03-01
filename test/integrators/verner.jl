using DelayDiffEq, DDEProblemLibrary
using Test

# simple problems
@testset "simple problems" begin
    prob_ip = prob_dde_constant_1delay_ip
    prob_scalar = prob_dde_constant_1delay_scalar
    ts = 0:0.1:10

    # Vern6
    println("Vern6")
    @testset "Vern6" begin
        alg = MethodOfSteps(Vern6())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 1.5e-3
        @test sol_ip.errors[:final] < 7.0e-6
        @test sol_ip.errors[:l2] < 1.0e-3

        sol_scalar = solve(prob_scalar, alg)

        # fails due to floating point issues
        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts) atol = 1.0e-5
    end

    # Vern7
    println("Vern7")
    @testset "Vern7" begin
        alg = MethodOfSteps(Vern7())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 6.0e-4
        @test sol_ip.errors[:final] < 3.5e-7
        @test sol_ip.errors[:l2] < 3.0e-4

        sol_scalar = solve(prob_scalar, alg)

        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts) atol = 1.0e-5
    end

    # Vern8
    println("Vern8")
    @testset "Vern8" begin
        alg = MethodOfSteps(Vern8())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 2.0e-3
        @test sol_ip.errors[:final] < 1.8e-5
        @test sol_ip.errors[:l2] < 8.8e-4

        sol_scalar = solve(prob_scalar, alg)

        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts) atol = 1.0e-5
    end

    # Vern9
    println("Vern9")
    @testset "Vern9" begin
        alg = MethodOfSteps(Vern9())
        sol_ip = solve(prob_ip, alg)

        @test sol_ip.errors[:l∞] < 1.5e-3
        @test sol_ip.errors[:final] < 3.8e-6
        @test sol_ip.errors[:l2] < 6.5e-4

        sol_scalar = solve(prob_scalar, alg)

        @test sol_ip(ts, idxs = 1) ≈ sol_scalar(ts) atol = 1.0e-5
    end
end

# model of Mackey and Glass
println("Mackey and Glass")
@testset "Mackey and Glass" begin
    prob = prob_dde_DDETST_A1

    # Vern6
    solve(prob, MethodOfSteps(Vern6()))

    # Vern7
    solve(prob, MethodOfSteps(Vern7()))

    # Vern8
    solve(prob, MethodOfSteps(Vern8()))

    # Vern9
    solve(prob, MethodOfSteps(Vern9()))
end
