using Test, DelayDiffEq, OrdinaryDiffEqSDIRK
using OrdinaryDiffEqCore: IController
using SciMLBase: ReturnCode

@testset "history updates refresh implicit solver caches" begin
    τ = 0.01

    function robertson!(du, u, h, p, t)
        delayed_u2 = h(p, t - τ; idxs = 2)
        reaction = 0.04 * u[1] - 10_000 * delayed_u2 * u[3]
        decay = 3.0e7 * u[2]^2
        du[1] = -reaction
        du[2] = reaction - decay
        du[3] = decay
        return nothing
    end

    history(p, t; idxs = nothing) = idxs === nothing ? [1.0, 0.0, 0.0] : 0.0
    prob = DDEProblem(robertson!, history, (0.0, 0.5); constant_lags = [τ])
    method = SDIRK2()
    controller = IController(method; qmax = 10, qmax_first_step = 10)
    sol = solve(
        prob, MethodOfSteps(method); controller, reltol = 0.1, abstol = 1.0e-4,
        dt = 1.0e-6, maxiters = 1000, save_everystep = false,
        timeseries_errors = false, dense_errors = false
    )

    @test sol.retcode == ReturnCode.Success
end
