using OrdinaryDiffEqMultirate, DiffEqDevTools, Test, LinearAlgebra

@testset "MRI-GARK ERK22 family" begin
    @testset "MRIGARKERK22a" begin
        @testset "Construction" begin
            @test MRIGARKERK22a(m = 10) == MRIGARKERK22a(10)
            @test MRIGARKERK22a(m = 20) == MRIGARKERK22a(20)
        end

        @testset "Scalar out-of-place" begin
            prob = SplitODEProblem(
                (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
            )
            sol = solve(prob, MRIGARKERK22a(m = 10), dt = 0.05, adaptive = false)
            @test abs(sol.u[end] - exp(-1.0)) < 1.0e-3

            sol_a = solve(prob, MRIGARKERK22a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
            @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
        end

        @testset "Vector in-place" begin
            f1!(du, u, p, t) = (du .= -0.9 .* u)
            f2!(du, u, p, t) = (du .= -0.1 .* u)
            u0 = [1.0, 2.0, 3.0]
            prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

            sol = solve(prob, MRIGARKERK22a(m = 10), dt = 0.05, adaptive = false)
            @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-2

            sol_a = solve(prob, MRIGARKERK22a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
            @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
        end

        @testset "Convergence" begin
            analytic(u0, p, t) = u0 * exp(-3t)
            prob = SplitODEProblem(
                SplitFunction(
                    (u, p, t) -> -2u, (u, p, t) -> -u; analytic = analytic
                ),
                1.0, (0.0, 1.0)
            )
            dts = 1 ./ 2 .^ (6:-1:2)
            sim = test_convergence(dts, prob, MRIGARKERK22a(m = 8))
            @test sim.𝒪est[:l∞] ≈ 2 atol = 0.3
        end
    end

    @testset "MRIGARKERK22b" begin
        @testset "Construction" begin
            @test MRIGARKERK22b(m = 10) == MRIGARKERK22b(10)
            @test MRIGARKERK22b(m = 20) == MRIGARKERK22b(20)
        end

        @testset "Scalar out-of-place" begin
            prob = SplitODEProblem(
                (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
            )
            sol = solve(prob, MRIGARKERK22b(m = 10), dt = 0.05, adaptive = false)
            @test abs(sol.u[end] - exp(-1.0)) < 1.0e-3

            sol_a = solve(prob, MRIGARKERK22b(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
            @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
        end

        @testset "Vector in-place" begin
            f1!(du, u, p, t) = (du .= -0.9 .* u)
            f2!(du, u, p, t) = (du .= -0.1 .* u)
            u0 = [1.0, 2.0, 3.0]
            prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

            sol = solve(prob, MRIGARKERK22b(m = 10), dt = 0.05, adaptive = false)
            @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-2

            sol_a = solve(prob, MRIGARKERK22b(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
            @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
        end

        @testset "Convergence" begin
            analytic(u0, p, t) = u0 * exp(-3t)
            prob = SplitODEProblem(
                SplitFunction(
                    (u, p, t) -> -2u, (u, p, t) -> -u; analytic = analytic
                ),
                1.0, (0.0, 1.0)
            )
            dts = 1 ./ 2 .^ (6:-1:2)
            sim = test_convergence(dts, prob, MRIGARKERK22b(m = 8))
            @test sim.𝒪est[:l∞] ≈ 2 atol = 0.3
        end
    end

end

@testset "MRIGARKERK33a" begin
    @testset "Construction" begin
        @test MRIGARKERK33a(m = 10) == MRIGARKERK33a(10)
        @test MRIGARKERK33a(m = 20) == MRIGARKERK33a(20)
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKERK33a(m = 10), dt = 0.05, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-4

        sol_a = solve(prob, MRIGARKERK33a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MRIGARKERK33a(m = 10), dt = 0.05, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-3

        sol_a = solve(prob, MRIGARKERK33a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
    end

    @testset "Convergence" begin
        analytic(u0, p, t) = u0 * exp(-3t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -2u, (u, p, t) -> -u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        dts = 1 ./ 2 .^ (6:-1:2)
        sim = test_convergence(dts, prob, MRIGARKERK33a(m = 8))
        @test sim.𝒪est[:l∞] ≈ 3 atol = 0.3
    end
end
