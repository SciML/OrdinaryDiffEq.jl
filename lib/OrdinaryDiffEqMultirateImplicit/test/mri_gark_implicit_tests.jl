using OrdinaryDiffEqMultirateImplicit, DiffEqDevTools, Test, LinearAlgebra

@testset "MRIGARKIRK21a" begin
    @testset "Construction" begin
        @test MRIGARKIRK21a(m = 10).m == 10
        @test MRIGARKIRK21a(m = 20).m == 20
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKIRK21a(m = 10), dt = 0.05, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-3

        sol_a = solve(prob, MRIGARKIRK21a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MRIGARKIRK21a(m = 10), dt = 0.05, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-2

        sol_a = solve(prob, MRIGARKIRK21a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
    end

    @testset "Stiff slow component" begin
        analytic(u0, p, t) = u0 * exp(-51t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -u, (u, p, t) -> -50u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKIRK21a(m = 4), dt = 0.02, adaptive = false)
        @test abs(sol.u[end] - exp(-51.0)) < 1.0e-3
        @test solve(prob, MRIGARKIRK21a(m = 4), reltol = 1.0e-6, abstol = 1.0e-8).retcode ==
            ReturnCode.Success
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
        sim = test_convergence(dts, prob, MRIGARKIRK21a(m = 8))
        @test sim.𝒪est[:l∞] ≈ 2 atol = 0.3
    end
end

@testset "MRIGARKESDIRK34a" begin
    @testset "Construction" begin
        @test MRIGARKESDIRK34a(m = 10).m == 10
        @test MRIGARKESDIRK34a(m = 20).m == 20
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKESDIRK34a(m = 10), dt = 0.05, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-3

        sol_a = solve(prob, MRIGARKESDIRK34a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MRIGARKESDIRK34a(m = 10), dt = 0.05, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-2

        sol_a = solve(prob, MRIGARKESDIRK34a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
    end

    @testset "Stiff slow component" begin
        analytic(u0, p, t) = u0 * exp(-51t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -u, (u, p, t) -> -50u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKESDIRK34a(m = 4), dt = 0.02, adaptive = false)
        @test abs(sol.u[end] - exp(-51.0)) < 1.0e-3
        @test solve(prob, MRIGARKESDIRK34a(m = 4), reltol = 1.0e-6, abstol = 1.0e-8).retcode ==
            ReturnCode.Success
    end

    @testset "Convergence" begin
        analytic(u0, p, t) = u0 * exp(-3t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -2u, (u, p, t) -> -u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        dts = 1 ./ 2 .^ (7:-1:3)
        sim = test_convergence(dts, prob, MRIGARKESDIRK34a(m = 8))
        @test sim.𝒪est[:l∞] ≈ 3 atol = 0.3
    end
end

@testset "MRIGARKESDIRK46a" begin
    @testset "Construction" begin
        @test MRIGARKESDIRK46a(m = 10).m == 10
        @test MRIGARKESDIRK46a(m = 20).m == 20
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKESDIRK46a(m = 10), dt = 0.05, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-3

        sol_a = solve(prob, MRIGARKESDIRK46a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MRIGARKESDIRK46a(m = 10), dt = 0.05, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-2

        sol_a = solve(prob, MRIGARKESDIRK46a(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
    end

    @testset "Stiff slow component" begin
        analytic(u0, p, t) = u0 * exp(-51t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -u, (u, p, t) -> -50u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        sol = solve(prob, MRIGARKESDIRK46a(m = 4), dt = 0.02, adaptive = false)
        @test abs(sol.u[end] - exp(-51.0)) < 1.0e-3
        @test solve(prob, MRIGARKESDIRK46a(m = 4), reltol = 1.0e-6, abstol = 1.0e-8).retcode ==
            ReturnCode.Success
    end

    @testset "Embedded error estimate" begin
        # The final coupling row is zero, so `u - z[s]` degenerates to 0 here: the
        # step size can only respond to the tolerance via the tableau embedding.
        f1!(du, u, p, t) = (du .= -2 .* u)
        f2!(du, u, p, t) = (du .= -u)
        prob = SplitODEProblem(f1!, f2!, [1.0], (0.0, 1.0))
        loose = solve(prob, MRIGARKESDIRK46a(m = 8), reltol = 1.0e-4, abstol = 1.0e-4)
        tight = solve(prob, MRIGARKESDIRK46a(m = 8), reltol = 1.0e-8, abstol = 1.0e-8)
        @test tight.stats.naccept > loose.stats.naccept
        @test abs(tight.u[end][1] - exp(-3.0)) < abs(loose.u[end][1] - exp(-3.0))
    end

    @testset "Convergence" begin
        analytic(u0, p, t) = u0 * exp(-3t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -2u, (u, p, t) -> -u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        dts = 1 ./ 2 .^ (8:-1:4)
        sim = test_convergence(dts, prob, MRIGARKESDIRK46a(m = 8))
        @test sim.𝒪est[:l∞] ≈ 4 atol = 0.3
    end
end
