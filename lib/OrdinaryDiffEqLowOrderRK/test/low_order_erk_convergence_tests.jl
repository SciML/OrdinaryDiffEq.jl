# This definitely needs cleaning
using OrdinaryDiffEqLowOrderRK, ODEProblemLibrary, DiffEqDevTools
using Test, Random
Random.seed!(100)

## Convergence Testing
dts1 = 1 .// 2 .^ (9:-1:5)
dts2 = 1 .// 2 .^ (7:-1:3)
dts3 = 1 .// 2 .^ (12:-1:7)
dts4 = 1 .// 2 .^ (5:-1:3)
dts5 = 1 .// 2 .^ (3:-1:1)
dts6 = 1 .// 10 .^ (5:-1:1)
testTol = 0.2

f = (u, p, t) -> sin(u)
prob_ode_nonlinear = ODEProblem(
    ODEFunction(
        f;
        analytic = (u0, p, t) -> 2 * acot(
            exp(-t) *
                cot(0.5)
        )
    ), 1.0,
    (0.0, 0.5)
)

@testset "Explicit Solver Convergence Tests ($(["out-of-place", "in-place"][i]))" for i in 1:2
    prob = (
        ODEProblemLibrary.prob_ode_linear,
        ODEProblemLibrary.prob_ode_2Dlinear,
    )[i]
    dts = 1 .// 2 .^ (8:-1:4)
    @info "Very low order"
    sim = test_convergence(dts, prob, Euler())
    @test sim.𝒪est[:final] ≈ 1 atol = testTol
    sim2 = test_convergence(dts, prob, Heun())
    @test sim2.𝒪est[:l∞] ≈ 2 atol = testTol
    sim2 = test_convergence(dts, prob, Ralston())
    @test sim2.𝒪est[:l∞] ≈ 2 atol = testTol
    sim2 = test_convergence(dts, prob, Midpoint())
    @test sim2.𝒪est[:l∞] ≈ 2 atol = testTol
    sim3 = test_convergence(dts, prob, RK4())
    @test sim3.𝒪est[:l∞] ≈ 4 atol = testTol

    sim3 = test_convergence(dts2, prob, RKO65())
    @test sim3.𝒪est[:l∞] ≈ 5 atol = testTol

    sim3 = test_convergence(dts4, prob, FRK65())
    @test sim3.𝒪est[:l∞] ≈ 6 atol = 0.6

    sim3 = test_convergence(dts, prob, RKM())
    @test sim3.𝒪est[:l∞] ≈ 4 atol = 0.2

    sim_ps6 = test_convergence(dts2, prob_ode_nonlinear, PSRK4p7q6())
    @test sim_ps6.𝒪est[:l∞] ≈ 4 atol = testTol

    sim_ps5 = test_convergence(dts2, prob_ode_nonlinear, PSRK3p6q5())
    @test sim_ps5.𝒪est[:l∞] ≈ 3 atol = testTol

    sim_ps4 = test_convergence(dts2, prob_ode_nonlinear, PSRK3p5q4())
    @test sim_ps4.𝒪est[:l∞] ≈ 3 atol = testTol

    sim_ms5 = test_convergence(dts2, prob, MSRK5())
    @test sim_ms5.𝒪est[:l∞] ≈ 5 atol = testTol

    sim_ms6 = test_convergence(dts4, prob, MSRK6())
    @test sim_ms6.𝒪est[:l∞] ≈ 6 atol = testTol

    sim_ms54 = test_convergence(dts2, prob, Stepanov5())
    @test sim_ms54.𝒪est[:l∞] ≈ 5 atol = 0.5

    sim4 = test_convergence(dts, prob, BS3())
    @test sim4.𝒪est[:l2] ≈ 3 atol = testTol

    sim4 = test_convergence(dts2, prob, SIR54())
    @test sim4.𝒪est[:l2] ≈ 4.4 atol = testTol

    sim2 = test_convergence(dts, prob, Alshina2())
    @test sim2.𝒪est[:l∞] ≈ 2 atol = testTol

    sim3 = test_convergence(dts, prob, Alshina3())
    @test sim3.𝒪est[:l∞] ≈ 3 atol = testTol

    sim6 = test_convergence(dts4, prob, Alshina6())
    @test sim6.𝒪est[:l∞] ≈ 6 atol = testTol

    sim160 = test_convergence(dts, prob, Anas5(w = 2))
    @test sim160.𝒪est[:l2] ≈ 4 atol = 2 * testTol
end
