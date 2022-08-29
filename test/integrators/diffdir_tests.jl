using Test, OrdinaryDiffEq

for Alg in (Rosenbrock23, TRBDF2)
    tspan = (0, 10.0)
    u0 = ones(2)
    tstops = range(0, 10, length = 1000)
    sol = solve(ODEProblem((du, u, p, t) -> du .= sin.(u), u0, tspan), Tsit5())
    @test_nowarn solve(ODEProblem((du, u, p, t) -> sol(du, t), u0, tspan), Alg(),
                       tstops = tstops)
    @test_nowarn solve(ODEProblem((du, u, p, t) -> sol(du, t), u0, (tspan[2], tspan[1])),
                       Alg(), tstops = tstops)

    for u0 in (1.0, ones(2))
        u0 = 1.0
        sol = solve(ODEProblem((u, p, t) -> sin.(u), u0, tspan), Tsit5())
        @test_nowarn solve(ODEProblem((u, p, t) -> sol(t), u0, tspan), Alg(),
                           tstops = tstops)
        @test_nowarn solve(ODEProblem((u, p, t) -> sol(t), u0, (tspan[2], tspan[1])), Alg(),
                           tstops = tstops)
    end
end
