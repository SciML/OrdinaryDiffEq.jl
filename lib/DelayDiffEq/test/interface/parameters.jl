using DelayDiffEq
using Test

# Test parameterized delayed logistic equation

# delayed logistic equation
f_inplace(du, u, h, p, t) = (du[1] = p[1] * u[1] * (1 - h(p, t - 1; idxs = 1)))
f_scalar(u, h, p, t) = p[1] * u * (1 - h(p, t - 1))

# simple history function
h(p, t; idxs = nothing) = 0.1

@testset for inplace in (true, false)
    # define problem
    # we specify order_discontinuity_t0 = 1 to indicate that the discontinuity at
    # t = 0 is of first order
    prob = DDEProblem(
        inplace ? f_inplace : f_scalar,
        inplace ? [0.1] : 0.1,
        h, (0.0, 50.0), [0.3];
        constant_lags = [1],
        order_discontinuity_t0 = 1
    )

    # solve problem with initial parameter:
    sol1 = solve(prob, MethodOfSteps(Tsit5()))
    @test length(sol1) == 21
    @test first(sol1(12)) ≈ 0.884 atol = 1.0e-4
    @test first(sol1.u[end]) ≈ 1 atol = 1.0e-5

    # solve problem with updated parameter
    prob.p[1] = 1.4
    sol2 = solve(prob, MethodOfSteps(Tsit5()))
    @test length(sol2) == 47
    @test first(sol2(12)) ≈ 1.125 atol = 5.0e-4
    @test first(sol2.u[end]) ≈ 0.994 atol = 2.0e-5
end
