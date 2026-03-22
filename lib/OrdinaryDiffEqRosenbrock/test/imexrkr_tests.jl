using OrdinaryDiffEqRosenbrock, DiffEqDevTools, SciMLBase, Test, LinearAlgebra

# Simple IMEX test problem: du/dt = -u (implicit) + 0 (explicit)
# Analytic solution: u(t) = u0 * exp(-t)
let
    f1_ip = (du, u, p, t) -> (du .= -u)
    f2_ip = (du, u, p, t) -> (du .= zero(u))
    f_ip = SplitFunction(f1_ip, f2_ip; analytic = (u0, p, t) -> u0 .* exp(-t))
    prob_ip = SplitODEProblem(f_ip, [1.0], (0.0, 1.0))

    f1_oop = (u, p, t) -> -u
    f2_oop = (u, p, t) -> zero(u)
    f_oop = SplitFunction(f1_oop, f2_oop; analytic = (u0, p, t) -> u0 .* exp(-t))
    prob_oop = SplitODEProblem(f_oop, [1.0], (0.0, 1.0))

    dts = (1 / 2) .^ (7:-1:4)
    testTol = 0.3

    @testset "IMEXRKR_3_2 convergence (in-place)" begin
        sim = test_convergence(dts, prob_ip, IMEXRKR_3_2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
    end

    @testset "IMEXRKR_3_2 convergence (out-of-place)" begin
        sim = test_convergence(dts, prob_oop, IMEXRKR_3_2())
        @test sim.𝒪est[:final] ≈ 2 atol = testTol
    end

    @testset "IMEXRKR_3_2 dense output" begin
        sol = solve(prob_ip, IMEXRKR_3_2(), dt = 0.01, adaptive = false, dense = true)
        @test SciMLBase.successful_retcode(sol)
        # Order-2 integrator + linear dense: global error O(dt²) ≈ 1e-4 at step boundaries;
        # linear interpolation within each step adds another O(dt²) ≈ 1e-4.
        for t_test in 0.005:0.1:0.995
            @test abs(sol(t_test)[1] - exp(-t_test)) < 1.0e-3
        end
    end
end
