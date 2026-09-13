using DelayDiffEq, OrdinaryDiffEqBDF, Test

@testset "NordsieckBDF trial-step history" begin
    # Exact solution exp(-t), including the history. DelayDiffEq copies the
    # trial polynomial before the acceptance controller runs, so publishing k
    # only on acceptance feeds the preceding step's polynomial into the RHS.
    lag = 0.2
    f(u, h, p, t) = -exp(-lag) * h(p, t - lag)
    function f!(du, u, h, p, t)
        du[1] = -exp(-lag) * h(p, t - lag; idxs = 1)
        return nothing
    end
    h(p, t) = exp(-t)
    hvec(p, t; idxs = nothing) = idxs === nothing ? [exp(-t)] : exp(-t)
    scalar(u) = u isa Number ? u : u[1]
    for iip in (false, true), with_callback in (false, true)

        @testset "inplace=$iip callback=$with_callback" begin
            prob = iip ?
                   DDEProblem(f!, [1.0], hvec, (0.0, 2.0); constant_lags = [lag]) :
                   DDEProblem(f, 1.0, h, (0.0, 2.0); constant_lags = [lag])
            events = Float64[]
            callback = with_callback ?
                       ContinuousCallback(
                (u, t, integrator) -> scalar(u) - 0.5,
                integrator -> push!(events, integrator.t)
            ) : nothing
            integrator = init(
                prob, MethodOfSteps(NordsieckBDF()); callback,
                dt = 0.08, dtmax = 0.08, abstol = 1.0e-10, reltol = 1.0e-9
            )
            solve!(integrator)
            sol = integrator.sol
            @test sol.retcode == ReturnCode.Success
            @test sol.stats.nreject > 0
            ts = range(0.01, 1.99; length = 201)
            @test maximum(t -> abs(scalar(sol(t)) - exp(-t)), ts) < 1.0e-6
            # Check the distinct history solution too: correct final-output
            # interpolation alone does not ensure correct delayed RHS values.
            history_sol = integrator.integrator.sol
            @test maximum(t -> abs(scalar(history_sol(t)) - exp(-t)), ts) < 1.0e-6
            if with_callback
                @test length(events) == 1
                @test isapprox(only(events), log(2.0); atol = 1.0e-6)
            end
        end
    end
end
