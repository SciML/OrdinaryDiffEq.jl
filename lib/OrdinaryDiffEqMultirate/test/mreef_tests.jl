using OrdinaryDiffEqMultirate, OrdinaryDiffEqLowOrderRK, DiffEqDevTools, Test, LinearAlgebra

@testset "MREEF" begin
    @testset "Construction" begin
        @test MREEF() == MREEF(4, 4, :harmonic)
        @test MREEF(m = 8, order = 3, seq = :romberg) == MREEF(8, 3, :romberg)
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )

        sol = solve(prob, MREEF(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-6

        sol_a = solve(prob, MREEF(m = 10, order = 4), reltol = 1.0e-8, abstol = 1.0e-8)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-6
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MREEF(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-6

        sol_a = solve(prob, MREEF(m = 10, order = 4), reltol = 1.0e-8, abstol = 1.0e-8)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-6
    end

    @testset "Order convergence" begin
        f1_analytic = (u0, p, t) -> u0 * exp(-t)
        f1_inplace = (du, u, p, t) -> (du .= -0.9 .* u)
        f2_inplace = (du, u, p, t) -> (du .= -0.1 .* u)

        f1_func = ODEFunction(f1_inplace; analytic = f1_analytic)
        prob = SplitODEProblem(f1_func, f2_inplace, [1.0], (0.0, 1.0))

        dts = 1 .// 2 .^ (6:-1:2)
        testTol = 0.3

        for target_order in [2, 3, 4]
            sim = test_convergence(
                dts, prob, MREEF(m = 10, order = target_order), adaptive = false
            )
            @test sim.𝒪est[:final] ≈ target_order atol = testTol
        end

        for target_order in [2, 3]
            sim = test_convergence(
                dts, prob, MREEF(m = 10, order = target_order, seq = :romberg),
                adaptive = false
            )
            @test sim.𝒪est[:final] ≈ target_order atol = testTol
        end
    end

    @testset "Nonlinear coupled system" begin
        function ff!(du, u, p, t)
            du[1] = 0.0
            du[2] = 20.0 * (2.0 * u[1] - u[1]^2 * u[2])
        end
        function fs!(du, u, p, t)
            du[1] = 1.0 + u[1]^2 * u[2] - 3.0 * u[1]
            du[2] = 0.0
        end
        prob = SplitODEProblem(ff!, fs!, [1.5, 3.0], (0.0, 0.5))
        ref = solve(prob, SplitEuler(), dt = 1.0e-6, adaptive = false)

        sol = solve(prob, MREEF(m = 20, order = 4), dt = 0.01, adaptive = false)
        @test norm(sol.u[end] - ref.u[end]) < 1.0e-3

        sol_a = solve(prob, MREEF(m = 20, order = 4), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - ref.u[end]) < 1.0e-3
    end

    @testset "Stats tracking" begin
        # For MREEF(m, order, seq=:harmonic): ns[j] = j
        # f1 evals per step = m * sum(1:order)
        # f2 evals per step = sum(1:order)
        # So nf / nf2 ≈ m (the extra init eval shifts it slightly)
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        prob = SplitODEProblem(f1!, f2!, [1.0], (0.0, 0.5))
        m_val = 5
        sol = solve(prob, MREEF(m = m_val, order = 3), dt = 0.1, adaptive = false)

        @test sol.stats.nf > 0
        @test sol.stats.nf2 > 0
        @test sol.stats.nf > sol.stats.nf2
        @test sol.stats.nf ≈ m_val * sol.stats.nf2 atol = m_val
    end

    @testset "Complex numbers" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9im .* u,
            (u, p, t) -> -0.1im .* u,
            [1.0 + 0.0im], (0.0, 1.0)
        )
        sol = solve(prob, MREEF(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test abs(sol.u[end][1] - exp(-1.0im)) < 1.0e-6
    end

    @testset "SplitFunction wrapper" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        ff = SplitFunction(f1!, f2!)
        prob = ODEProblem(ff, [1.0, 2.0], (0.0, 1.0))
        sol = solve(prob, MREEF(m = 10, order = 4), dt = 0.1, adaptive = false)
        @test norm(sol.u[end] - [1.0, 2.0] .* exp(-1.0)) < 1.0e-6
    end
end
