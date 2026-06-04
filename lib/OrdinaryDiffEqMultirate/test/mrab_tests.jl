using OrdinaryDiffEqMultirate, DiffEqDevTools, Test, LinearAlgebra

@testset "MRAB" begin
    @testset "Construction" begin
        @test MRAB() == MRAB(2, 4)
        @test MRAB(k = 3, m = 8) == MRAB(3, 8)
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )
        for k in 1:4
            sol = solve(prob, MRAB(k = k, m = 10), dt = 0.02, adaptive = false)
            @test abs(sol.u[end] - exp(-1.0)) < 2.0e-2
        end

        sol_a = solve(prob, MRAB(k = 2, m = 10), reltol = 1.0e-8, abstol = 1.0e-8)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        for k in 1:4
            sol = solve(prob, MRAB(k = k, m = 10), dt = 0.02, adaptive = false)
            @test norm(sol.u[end] - u0 .* exp(-1.0)) < 5.0e-2
        end

        sol_a = solve(prob, MRAB(k = 2, m = 10), reltol = 1.0e-8, abstol = 1.0e-8)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
    end

    @testset "Convergence" begin
        # Order is globally 1 (frozen-slow Lie splitting); higher k improves the
        # error constant but not the asymptotic rate.
        analytic(u0, p, t) = u0 * exp(-t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        dts = 1 ./ 2 .^ (8:-1:4)
        for k in 1:4
            sim = test_convergence(dts, prob, MRAB(k = k, m = 4))
            @test sim.𝒪est[:l∞] ≈ 1 atol = 0.2
        end
    end
end
