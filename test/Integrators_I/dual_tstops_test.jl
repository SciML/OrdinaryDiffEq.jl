using OrdinaryDiffEq, Test, ForwardDiff

# Dual tstops with a Float64 tspan used to throw a MethodError
# https://github.com/SciML/OrdinaryDiffEq.jl/issues/4989
@testset "Dual tstops with Float64 tspan" begin
    prob = ODEProblem((u, p, t) -> -p[1] * u, 1.0, (0.0, 2.0), [0.5])

    f(ts) = solve(prob, Tsit5(); tstops = [ts])(2.0)
    # reference: promote tspan to the Dual type by hand
    g(ts) = solve(remake(prob; tspan = (zero(ts), 2one(ts))), Tsit5(); tstops = [ts])(2.0)

    @test ForwardDiff.derivative(f, 1.0) ≈ ForwardDiff.derivative(g, 1.0)
    # the exact solution doesn't depend on tstops, so the derivative is ~0
    @test abs(ForwardDiff.derivative(f, 1.0)) < 1.0e-5
end

@testset "Dual tstops with a stopping callback" begin
    u0 = 1.0
    function h(tstop)
        cb = DiscreteCallback(
            (u, t, integrator) -> t == tstop,
            integrator -> (integrator.p = zero(integrator.p))
        )
        prob = ODEProblem((u, p, t) -> p * u, u0, (0.0, 1.0), 1.0)
        return solve(prob, Tsit5(); callback = cb, tstops = [tstop])(1.0)
    end
    # u' = u until tstop, then u' = 0, so u(1) = exp(tstop) * u0
    @test ForwardDiff.derivative(h, 0.5) ≈ exp(0.5) * u0
end
