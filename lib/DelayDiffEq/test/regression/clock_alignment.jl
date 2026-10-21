using DelayDiffEq, OrdinaryDiffEqTsit5, Test
using OrdinaryDiffEqCore: PIController
using SciMLBase: ReturnCode

@testset "DDE clock and state intervals" begin
    history(p, t) = p
    f(u, h, p, t) = h(p, t - 100.0)
    for t0 in (0.0, 2.0^40)
        integrators = map((1.0, -2.0)) do slope
            prob = DDEProblem(
                f, 1.0, history, (t0, t0 + 10.0), slope; constant_lags = [100.0]
            )
            init(
                prob, MethodOfSteps(Tsit5()); dt = 0.1,
                controller = PIController(
                    0.0, 0.0; qmin = 1.0, qmax = 1.0, qmax_first_step = 1.0
                ),
                abstol = 1.0e-10, reltol = 1.0e-10
            )
        end
        foreach(step!, integrators)
        first_steps = map(integrator -> (integrator.t, integrator.u), integrators)
        for i in (2, 1)
            integrator = integrators[i]
            sol = solve!(integrator)
            t1, u1 = first_steps[i]
            slope = sol.prob.p
            @test sol.retcode == ReturnCode.Success
            @test sol.t[end] == t0 + 10.0
            # The initial dt is chosen before calc_dt_propose! aligns subsequent steps.
            @test maximum(
                abs(u - u1 - slope * (t - t1))
                    for (t, u) in zip(sol.t[2:end], sol.u[2:end])
            ) <= 1.0e-9
            if i == 2
                @test (integrators[1].t, integrators[1].u) == first_steps[1]
            end
        end
    end
end

@testset "DDE interpolated history" begin
    history(p, t) = 1.0
    f(u, h, p, t) = h(p, t - 1.0)
    prob = DDEProblem(f, 1.0, history, (0.0, 3.0); constant_lags = [1.0])
    reference(t) = 1 + t + max(t - 1, 0)^2 / 2 + max(t - 2, 0)^3 / 6
    sol = solve(prob, MethodOfSteps(Tsit5()); abstol = 1.0e-11, reltol = 1.0e-11)
    saved_t, saved_u = copy(sol.t), copy(sol.u)
    @test sol.retcode == ReturnCode.Success
    @test sol.t[end] == 3.0
    @test maximum(abs(sol(t) - reference(t)) for t in 0.0:0.125:3.0) < 1.0e-9
    repeated = solve(prob, MethodOfSteps(Tsit5()); abstol = 1.0e-11, reltol = 1.0e-11)
    @test repeated.t == saved_t
    @test repeated.u == saved_u
    @test sol.t == saved_t
    @test sol.u == saved_u
end
