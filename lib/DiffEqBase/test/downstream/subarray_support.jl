using OrdinaryDiffEq
using Test

# Test for SubArray support with ODEProblem
# Regression test for https://github.com/SciML/OrdinaryDiffEq.jl/issues/2900
# Fixed in https://github.com/SciML/DiffEqBase.jl/pull/1219

@testset "SubArray as initial conditions" begin
    # Test 1: Simple ODE with SubArray initial conditions
    @testset "Basic ODE with SubArray" begin
        f(du, u, p, t) = du .= u
        u0_full = ones(10)
        u0 = @view u0_full[1:5]

        # This should not throw NoFunctionWrapperFoundError
        prob = ODEProblem(f, u0, (0.0, 1.0))
        sol = solve(prob, Tsit5())

        @test sol.retcode == ReturnCode.Success
        @test length(sol.u[end]) == 5
    end

    # Test 2: Simple pendulum with SubArray initial conditions
    @testset "Simple pendulum with SubArray" begin
        g = 9.81
        L = 1.0

        # Initial Conditions as SubArray
        u₀ = @view [0, π / 60][:]  # Initial speed and initial angle
        tspan = (0.0, 6.3)

        # Define the pendulum problem
        function simplependulum(du, u, p, t)
            θ = u[1]
            dθ = u[2]
            du[1] = dθ
            du[2] = -(g / L) * θ
        end

        # This should not throw NoFunctionWrapperFoundError
        prob = ODEProblem(simplependulum, u₀, tspan)
        integrator = init(prob, Tsit5(); reltol = 1.0e-6)

        @test integrator !== nothing

        # Solve and verify the solution
        sol = solve(prob, Tsit5(); reltol = 1.0e-6)
        @test sol.retcode == ReturnCode.Success
        @test length(sol.u[end]) == 2
    end

    # Test 3: Various SubArray slicing patterns
    @testset "Different SubArray patterns" begin
        f(du, u, p, t) = du .= 2u

        # Test different ways to create SubArrays
        @testset "range slice" begin
            arr = [1.0, 2.0, 3.0]
            u0 = @view arr[1:2]
            prob = ODEProblem(f, u0, (0.0, 1.0))
            sol = solve(prob, Tsit5())
            @test sol.retcode == ReturnCode.Success
        end

        @testset "colon slice" begin
            arr = [1.0, 2.0, 3.0]
            u0 = @view arr[:]
            prob = ODEProblem(f, u0, (0.0, 1.0))
            sol = solve(prob, Tsit5())
            @test sol.retcode == ReturnCode.Success
        end

        @testset "offset range" begin
            arr = [1.0, 2.0, 3.0]
            u0 = @view arr[2:3]
            prob = ODEProblem(f, u0, (0.0, 1.0))
            sol = solve(prob, Tsit5())
            @test sol.retcode == ReturnCode.Success
        end

        @testset "view of constructed array" begin
            arr = ones(5)
            u0 = @view arr[1:3]
            prob = ODEProblem(f, u0, (0.0, 1.0))
            sol = solve(prob, Tsit5())
            @test sol.retcode == ReturnCode.Success
        end
    end
end
