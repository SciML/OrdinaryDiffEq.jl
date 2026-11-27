using DelayDiffEq, DDEProblemLibrary
using OrdinaryDiffEqLowOrderRK
using Logging, Test

const prob = prob_dde_constant_1delay_ip
const alg = MethodOfSteps(BS3(); constrained = false)

intermediate_fracs(logger) = [
    r.kwargs[:progress] for r in logger.logs
        if haskey(r.kwargs, :progress) && r.kwargs[:progress] isa Real &&
        !iszero(r.kwargs[:progress])
]

function progress_solve!(integrator)
    logger = Test.TestLogger(min_level = LogLevel(-1))
    integrator.opts.progress_min_interval_ns = 0  # throttle open: only the tspan gate decides
    with_logger(() -> solve!(integrator), logger)
    return intermediate_fracs(logger)
end

@testset "progress keywords are passed through" begin
    @test init(prob, alg; dt = 0.01).opts.progress_steps == 1
    @test init(prob, alg; dt = 0.01, progress_min_interval = 0.05).opts.progress_min_interval_ns ==
        50_000_000
    @test_throws ArgumentError init(prob, alg; dt = 0.01, progress_tspan = -1)
    integrator = init(prob, alg; dt = 0.01, progress = true, progress_tspan = 0.1)
    @test integrator.opts.progress_tspan == 0.1
    fracs = progress_solve!(integrator)
    @test !isempty(fracs)
    @test all(diff(fracs) .>= 0.1 * (1 - 1.0e-3))
    @test length(fracs) <= 10
end

@testset "reinit! resets the progress state" begin
    integrator = init(prob, alg; dt = 0.01, progress = true)
    first_fracs = progress_solve!(integrator)
    @test integrator.opts.progress_last_frac > 0.9
    last_ns = integrator.opts.progress_last_ns
    reinit!(integrator)
    @test integrator.opts.progress_last_frac == 0
    @test integrator.opts.progress_last_ns >= last_ns  # clock restarted
    @test length(progress_solve!(integrator)) == length(first_fracs)
end
