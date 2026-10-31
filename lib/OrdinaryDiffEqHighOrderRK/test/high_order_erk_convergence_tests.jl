# This definitely needs cleaning
using OrdinaryDiffEqHighOrderRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random, SciMLBase
Random.seed!(100)

dts5 = 1 .// 2 .^ (3:-1:1)
testTol = 0.2

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]

    sim3 = test_convergence(dts5, prob, PFRK87())
    @test sim3.𝒪est[:l∞] ≈ 8.4 atol = 0.2
end

@testset "PFRK87 adaptive error estimate (#4403)" begin
    prob_const = ODEProblem((u, p, t) -> one(u), 0.0, (0.0, 1.0))
    sol = solve(prob_const, PFRK87(); abstol = 1.0e-6, reltol = 1.0e-6)
    @test SciMLBase.successful_retcode(sol)
    @test sol.stats.naccept < 20
    @test sol.u[end] ≈ 1.0

    prob_decay = ODEProblem((u, p, t) -> -u, 1.0, (0.0, 1.0))
    errs = map((1.0e-4, 1.0e-8)) do tol
        s = solve(prob_decay, PFRK87(); abstol = tol, reltol = tol)
        @test SciMLBase.successful_retcode(s)
        abs(s.u[end] - exp(-1.0))
    end
    @test errs[2] < 1.0e-9
end
