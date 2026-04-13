using OrdinaryDiffEqSIMDRK, DiffEqDevTools, Test

function nonauto1(u, p, t)
    x, _ = u
    return [t * x, 0]
end

function nonauto2(u, p, t)
    _, y = u
    return [y, t * y]
end

function analytic(u0, p, t)
    x0, y0 = u0
    et = exp(t^2 / 2)
    return [et * (x0 + t * y0), et * y0]
end

u0 = [1.1, 2.2]
tspan = (0.0, 1.0)
prob1 = ODEProblem(
    ODEFunction{true}(
        (du, u, p, t) -> du .= nonauto1(u, p, t) .+
            nonauto2(u, p, t),
        analytic = analytic
    ),
    u0, tspan
)
prob2 = ODEProblem(
    ODEFunction{false}(
        (u, p, t) -> nonauto1(u, p, t) .+ nonauto2(u, p, t),
        analytic = analytic
    ),
    u0, tspan
)

for prob in [prob2]
    #=prob1,=#
    dts1 = 1 .// 2 .^ (7:-1:4)
    dts2 = 1 .// 2 .^ (4:-1:1)

    sim = test_convergence(dts1, prob, MER5v2())
    @test 5 <= sim.𝒪est[:l∞] <= 6

    sim = test_convergence(dts1, prob, MER6v2())
    @test 6 <= sim.𝒪est[:l∞] <= 7

    sim = test_convergence(dts2, prob, RK6v4())
    @test sim.𝒪est[:l∞] ≈ 6 atol = 1.0e-2
end
