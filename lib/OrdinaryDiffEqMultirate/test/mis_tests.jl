using OrdinaryDiffEqMultirate, DiffEqDevTools, Test, LinearAlgebra

@testset "MIS" begin
    @testset "Construction" begin
        @test MIS(m = 10) == MIS(10)
        @test MIS(m = 20) == MIS(20)
    end

    @testset "Scalar out-of-place" begin
        prob = SplitODEProblem(
            (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u, 1.0, (0.0, 1.0)
        )
        sol = solve(prob, MIS(m = 10), dt = 0.05, adaptive = false)
        @test abs(sol.u[end] - exp(-1.0)) < 1.0e-3

        sol_a = solve(prob, MIS(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test abs(sol_a.u[end] - exp(-1.0)) < 1.0e-3
    end

    @testset "Vector in-place" begin
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        u0 = [1.0, 2.0, 3.0]
        prob = SplitODEProblem(f1!, f2!, u0, (0.0, 1.0))

        sol = solve(prob, MIS(m = 10), dt = 0.05, adaptive = false)
        @test norm(sol.u[end] - u0 .* exp(-1.0)) < 1.0e-2

        sol_a = solve(prob, MIS(m = 10), reltol = 1.0e-6, abstol = 1.0e-6)
        @test norm(sol_a.u[end] - u0 .* exp(-1.0)) < 1.0e-3
    end

    @testset "Convergence" begin
        analytic(u0, p, t) = u0 * exp(-t)
        prob = SplitODEProblem(
            SplitFunction(
                (u, p, t) -> -0.9 * u, (u, p, t) -> -0.1 * u; analytic = analytic
            ),
            1.0, (0.0, 1.0)
        )
        dts = 1 ./ 2 .^ (6:-1:2)
        sim = test_convergence(dts, prob, MIS(m = 4))
        @test sim.𝒪est[:l∞] ≈ 2 atol = 0.3
    end
end
