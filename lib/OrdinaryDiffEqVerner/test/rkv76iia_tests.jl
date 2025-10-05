using OrdinaryDiffEqVerner, OrdinaryDiffEqCore, DiffEqBase, Test
using LinearAlgebra
using OrdinaryDiffEqSSPRK, DiffEqDevTools, Test, Random
import OrdinaryDiffEqLowStorageRK
import ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear

testTol=0.3

# Simple convergence test for RKV76IIa
@testset "RKV76IIa Convergence Tests" begin
    # Test problem: u' = cos(t), y(0) = 0
    # Exact solution: u(t) = sin(t)
    f = (u, p, t) -> cos(t)
    prob_ode_sin = ODEProblem(ODEFunction(f; analytic = (u0, p, t) -> sin(t)), 0.0, (0.0, 1.0))

    f = (du, u, p, t) -> du[1] = cos(t)
    prob_ode_sin_inplace = ODEProblem(ODEFunction(f; analytic = (u0, p, t) -> [sin(t)]), [0.0],
        (0.0, 1.0))

    # Test problem: u' = sin(u), y(0) = 1
    # Exact solution: u(t) = 2 * arccot(e^(-t) * cot(0.5))
    f = (u, p, t) -> sin(u)
    prob_ode_nonlinear = ODEProblem(
        ODEFunction(f;
            analytic = (u0, p, t) -> 2 * acot(exp(-t) *
                                            cot(0.5))), 1.0,
        (0.0, 0.5))

    f = (du, u, p, t) -> du[1] = sin(u[1])
    prob_ode_nonlinear_inplace = ODEProblem(
        ODEFunction(f;
            analytic = (u0, p, t) -> [
                2 * acot(exp(-t) * cot(0.5))
            ]),
        [1.0], (0.0, 0.5))

    test_problems_only_time = [prob_ode_sin, prob_ode_sin_inplace]
    # test_problems_linear = [prob_ode_linear, prob_ode_2Dlinear, prob_ode_bigfloat2Dlinear]
    # test_problems_nonlinear = [prob_ode_nonlinear, prob_ode_nonlinear_inplace]

    function f!(du, u, p, t)
        du[1] = -u[1]
    end

    function f(u, p, t)
        -u
    end

    # Out-of-place test
    prob_oop = ODEProblem(f, 1.0, (0.0, 1.0))
    sol_oop = solve(prob_oop, RKV76IIa(), abstol=1e-12, reltol=1e-12)
    @test sol_oop[end] ≈ exp(-1.0) atol=1e-10
    println("Out-of-place solution at t=1: ", sol_oop[end])
    println("Expected value: ", exp(-1.0))

    dts = 1 .// 2 .^ (2:6)  # [1/16, 1/32, 1/64, 1/128, 1/256]
    # Tried other ranges of dts
    # dts = (1/2) .^ (5:9)
    # dts = 1 .// 2 .^ (4:8)

    errors = zeros(length(dts))
    println("Testing order 7:")
    for (i, dt) in enumerate(dts)
        sol = solve(prob_oop, RKV76IIa(), dt=dt, adaptive=false)
        errors[i] = abs(sol[end] - exp(-1.0))
        println("dt = ", dt, ", error = ", errors[i])
    end

    # Check convergence order, this is failing for now
    for i in 2:length(errors)
        order = log(errors[i-1]/errors[i]) / log(2)
        println("Order between dt=", dts[i-1], " and dt=", dts[i], ": ", order)
        @test order ≈ 7 atol=testTol  # Allow some tolerance

    end

    # Additional tests disabled for now
    # println("RKV76IIa")
    # alg = RKV76IIa()
    # for prob in test_problems_only_time
    #     sim = test_convergence(dts, prob, alg)
    #     @test sim.𝒪est[:final]≈OrdinaryDiffEqVerner.alg_order(alg) atol=testTol
    # end
    # for prob in test_problems_linear
    #     sim = test_convergence(dts, prob, alg)
    #     @test sim.𝒪est[:final]≈OrdinaryDiffEqVerner.alg_order(alg) atol=testTol
    # end
    # for prob in test_problems_nonlinear
    #     sim = test_convergence(dts, prob, alg)
    #     @test sim.𝒪est[:final]≈OrdinaryDiffEqVerner.alg_order(alg) atol=testTol
    # end
end

# # Test with a system of ODEs
# @testset "RKV76IIa System Test" begin
#     # Lorenz system
#     function lorenz!(du, u, p, t)
#         σ, ρ, β = p
#         du[1] = σ * (u[2] - u[1])
#         du[2] = u[1] * (ρ - u[3]) - u[2]
#         du[3] = u[1] * u[2] - β * u[3]
#     end

#     u0 = [1.0, 0.0, 0.0]
#     p = [10.0, 28.0, 8/3]
#     tspan = (0.0, 0.1)
#     prob = ODEProblem(lorenz!, u0, tspan, p)

#     sol = solve(prob, RKV76IIa(), abstol=1e-10, reltol=1e-10)
#     @test length(sol.t) > 2  # Should have taken multiple steps
#     @test all(isfinite.(sol[end]))  # Solution should be finite
# end

# # Test lazy vs non-lazy
# @testset "RKV76IIa Lazy Tests" begin
#     prob = ODEProblem((u,p,t) -> -u, 1.0, (0.0, 1.0))

#     sol_lazy = solve(prob, RKV76IIa(lazy=true), abstol=1e-10, reltol=1e-10)
#     sol_not_lazy = solve(prob, RKV76IIa(lazy=false), abstol=1e-10, reltol=1e-10)

#     # Both should give the same result
#     @test sol_lazy[end] ≈ sol_not_lazy[end] atol=1e-12
# end
# # In-place test
    # prob_ip = ODEProblem(f!, [1.0], (0.0, 1.0))
    # sol_ip = solve(prob_ip, RKV76IIa(), abstol=1e-12, reltol=1e-12)
    # @test sol_ip[end][1] ≈ exp(-1.0) atol=1e-10
    # println("In-place solution at t=1: ", sol_ip[end][1])
    # println("Expected value: ", exp(-1.0))

    # Test that it's order 7
    # dts = (1/2) .^ (9:5)
    # dts = 1 .// 2 .^ (8:-1:4)
