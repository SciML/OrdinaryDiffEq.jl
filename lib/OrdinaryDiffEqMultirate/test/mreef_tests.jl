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
        # For MREEF(m, order, seq=:harmonic): ns[j] = j, so the step itself costs
        # m * sum(1:order) evaluations of f1 and sum(1:order) of f2. On top of that
        # every step evaluates the right endpoint derivative once for the Hermite
        # interpolant, and initialize! evaluates the left endpoint once per solve.
        # (The left endpoint of later steps is reused from the first substep, free.)
        f1!(du, u, p, t) = (du .= -0.9 .* u)
        f2!(du, u, p, t) = (du .= -0.1 .* u)
        prob = SplitODEProblem(f1!, f2!, [1.0], (0.0, 0.5))
        m_val, order = 5, 3
        sol = solve(prob, MREEF(m = m_val, order = order), dt = 0.1, adaptive = false)
        nsteps = length(sol.t) - 1

        @test sol.stats.nf > sol.stats.nf2
        @test sol.stats.nf == (m_val * sum(1:order) + 1) * nsteps + 1
        @test sol.stats.nf2 == (sum(1:order) + 1) * nsteps + 1
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

# The Hermite interpolant is built from the endpoint derivatives k[1] and k[2].
# None of these methods are FSAL, so unless perform_step! writes both, sol(t)
# between step nodes silently degrades to the accuracy of a stale slope: the
# interpolation error used to sit at ~7e-3 no matter how small the node error was.
@testset "Dense output tracks node accuracy" begin
    w, e = 50.0, 0.1
    exact(t) = exp(-e * t) * [cos(w * t), sin(w * t)]
    f1!(du, u, p, t) = (du[1] = -w * u[2]; du[2] = w * u[1]; nothing)
    f2!(du, u, p, t) = (du .= -e .* u; nothing)
    f1(u, p, t) = [-w * u[2], w * u[1]]
    f2(u, p, t) = -e .* u

    probs = (
        SplitODEProblem(f1!, f2!, [1.0, 0.0], (0.0, 1.0)),
        SplitODEProblem(f1, f2, [1.0, 0.0], (0.0, 1.0)),
    )
    algs = (MREEF(m = 8, order = 4), MRAB(k = 2, m = 8), MRIGARKERK45a(m = 8), MIS(m = 8))

    for prob in probs, alg in algs
        sol = solve(prob, alg, dt = 1.0e-3, adaptive = false)
        node = maximum(
            norm(sol.u[i] - exact(sol.t[i])) / norm(exact(sol.t[i]))
                for i in 2:length(sol.t)
        )
        interp = maximum(
            begin
                    tm = (sol.t[i] + sol.t[i + 1]) / 2
                    norm(sol(tm) - exact(tm)) / norm(exact(tm))
                end for i in 1:(length(sol.t) - 1)
        )
        # Interpolating must not be orders of magnitude worse than landing on a node.
        @test interp < max(100 * node, 1.0e-7)
    end
end
