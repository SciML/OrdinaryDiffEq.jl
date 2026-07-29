using DelayDiffEq
using SciMLBase: DDEProblem, init, solve!, successful_retcode
using Test

@testset "DelayDiffEq 6.1.0 DEOptions layout" begin
    f(du, u, h, p, t) = (du[1] = p * u[1] * (1 - h(p, t - 1; idxs = 1)))
    h(p, t; idxs = nothing) = 0.1
    prob = DDEProblem(f, [0.1], h, (0.0, 2.0), 0.3; constant_lags = [1])

    @test Base.pkgversion(DelayDiffEq) == v"6.1.0"
    integrator = init(prob)
    @test hasproperty(integrator.opts, :stage_limiter!)
    @test hasproperty(integrator.opts, :step_limiter!)
    @test successful_retcode(solve!(integrator))
end
