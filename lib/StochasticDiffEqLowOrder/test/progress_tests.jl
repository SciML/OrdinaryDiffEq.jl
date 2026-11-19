using StochasticDiffEqLowOrder, Test
using Base.CoreLogging: LogLevel, with_logger  # Logging isn't a test dependency here

intermediate_fracs(logger) = [
    r.kwargs[:progress] for r in logger.logs
        if haskey(r.kwargs, :progress) && r.kwargs[:progress] isa Real &&
        !iszero(r.kwargs[:progress])
]

prob = SDEProblem((u, p, t) -> -u, (u, p, t) -> 0.1u, 1.0, (0.0, 10.0))

@testset "Defaults and keywords reach the SDE integrator" begin
    integrator = init(prob, EM(); dt = 1.0e-3)
    @test integrator.opts.progress_steps == 1
    integrator = init(
        prob, EM(); dt = 1.0e-3, progress = true, progress_tspan = 0.1, progress_min_interval = 0
    )
    @test integrator.opts.progress_tspan == 0.1
    @test integrator.opts.progress_min_interval_ns == 0
end

@testset "tspan gate on an SDE solve" begin
    logger = Test.TestLogger(min_level = LogLevel(-1))
    with_logger(logger) do
        solve(prob, EM(); dt = 1.0e-3, progress = true, progress_min_interval = 0)
    end
    fracs = intermediate_fracs(logger)
    @test 90 <= length(fracs) <= 100
    @test all(diff(fracs) .>= 0.01 * (1 - 1.0e-3))
    @test last(logger.logs).kwargs[:progress] == "done"
end
