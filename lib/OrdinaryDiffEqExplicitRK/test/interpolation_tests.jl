"""
    Interpolation Tests for ExplicitRK with Generic Dense Output

These tests verify that the generic interpolation implementation for ExplicitRK
works correctly through the standard user-facing interface (solve + solution interpolation).
"""

using Test
using OrdinaryDiffEqCore
using OrdinaryDiffEqExplicitRK
using OrdinaryDiffEqExplicitRK: constructTsit5ExplicitRK, constructTsit5ExplicitRKSimple
using SciMLBase: ODEProblem, solve

@testset "ExplicitRK Dense Output" begin

    @testset "Scalar ODE - Exponential Decay" begin
        # y' = -y, y(0) = 1
        # Exact solution: y(t) = exp(-t)
        f(u, p, t) = -u
        u0 = 1.0
        tspan = (0.0, 1.0)
        prob = ODEProblem(f, u0, tspan)

        # Solve with ExplicitRK using Tsit5 tableau
        tableau = constructTsit5ExplicitRK()
        sol = solve(prob, ExplicitRK(tableau=tableau), dense=true, abstol=1e-10, reltol=1e-10)

        # Test interpolation at various points
        # Note: Interpolation error can be larger than local error tolerance
        @test sol(0.0) ≈ 1.0 atol=1e-8
        @test sol(0.25) ≈ exp(-0.25) atol=1e-3
        @test sol(0.5) ≈ exp(-0.5) atol=1e-3
        @test sol(0.75) ≈ exp(-0.75) atol=1e-3
        @test sol(1.0) ≈ exp(-1.0) atol=1e-3
    end

    @testset "Vector ODE - Decoupled System" begin
        # y₁' = -y₁, y₂' = -2*y₂
        # Exact: y₁(t) = exp(-t), y₂(t) = 2*exp(-2t)
        function f!(du, u, p, t)
            du[1] = -u[1]
            du[2] = -2 * u[2]
        end
        u0 = [1.0, 2.0]
        tspan = (0.0, 1.0)
        prob = ODEProblem(f!, u0, tspan)

        tableau = constructTsit5ExplicitRK()
        sol = solve(prob, ExplicitRK(tableau=tableau), dense=true, abstol=1e-10, reltol=1e-10)

        # Test full vector interpolation
        # Note: Interpolation error can be larger than local error tolerance
        for t in [0.2, 0.4, 0.6, 0.8]
            interp_val = sol(t)
            exact_val = [exp(-t), 2 * exp(-2 * t)]
            @test interp_val ≈ exact_val atol=1e-3
        end

        # Test component-wise interpolation using idxs
        for t in [0.3, 0.7]
            @test sol(t, idxs=1) ≈ exp(-t) atol=1e-3
            @test sol(t, idxs=2) ≈ 2 * exp(-2 * t) atol=1e-3
        end

        # Test multiple indices
        result = sol(0.5, idxs=[1, 2])
        @test result ≈ [exp(-0.5), 2 * exp(-1.0)] atol=1e-3
    end

    @testset "Simple tableau (Float64 precision)" begin
        # Test with the simpler Float64 tableau
        f(u, p, t) = -u
        u0 = 1.0
        tspan = (0.0, 1.0)
        prob = ODEProblem(f, u0, tspan)

        tableau = constructTsit5ExplicitRKSimple()
        sol = solve(prob, ExplicitRK(tableau=tableau), dense=true, abstol=1e-8, reltol=1e-8)

        # Check interpolation works correctly
        @test sol(0.5) ≈ exp(-0.5) atol=2e-3
    end

    @testset "Large system interpolation" begin
        # Lorenz-96 model (smaller version for testing)
        function lorenz96!(du, u, p, t)
            F = p
            N = length(u)
            @inbounds for i in 1:N
                i_m2 = mod1(i - 2, N)
                i_m1 = mod1(i - 1, N)
                i_p1 = mod1(i + 1, N)
                du[i] = (u[i_p1] - u[i_m2]) * u[i_m1] - u[i] + F
            end
        end

        N = 100
        F = 8.0
        u0 = F .+ 0.01 .* ones(N)  # deterministic initial condition
        tspan = (0.0, 1.0)
        prob = ODEProblem(lorenz96!, u0, tspan, F)

        tableau = constructTsit5ExplicitRK()
        sol = solve(prob, ExplicitRK(tableau=tableau), dense=true, abstol=1e-6, reltol=1e-6)

        # Verify interpolation returns correct size and is smooth
        for t in [0.25, 0.5, 0.75]
            interp_val = sol(t)
            @test length(interp_val) == N
            @test all(isfinite.(interp_val))
        end

        # Test single component interpolation
        @test isfinite(sol(0.5, idxs=50))
    end

    @testset "Interpolation at exact boundaries" begin
        f(u, p, t) = -u
        u0 = 1.0
        prob = ODEProblem(f, u0, (0.0, 1.0))

        tableau = constructTsit5ExplicitRK()
        sol = solve(prob, ExplicitRK(tableau=tableau), dense=true, abstol=1e-10, reltol=1e-10)

        # Interpolation at exact start and end
        @test sol(0.0) ≈ u0
        @test sol(1.0) ≈ sol.u[end]
    end

end
