# AD through integrator steps that send progress updates. The progress gate reads the
# clock (`time_ns()`) and logs, so it must be opaque to every AD mode: Mooncake forward
# mode, for instance, has no rule for the clock's foreign call.
using OrdinaryDiffEqTsit5, ForwardDiff, Test
using Base.CoreLogging: LogLevel, with_logger  # Logging isn't a dependency of this env
using ADTypes
import DifferentiationInterface as DI
using Mooncake  # Load Mooncake after DI to ensure extension is loaded
import Enzyme

function progress_steps_solution(p)
    prob = ODEProblem((u, p, t) -> -p * u, 1.0, (0.0, 1.0), p)
    integ = init(
        prob, Tsit5(); progress = true, save_everystep = false, dt = 0.1, adaptive = false
    )
    integ.opts.progress_min_interval_ns = 0  # throttle open: every gate pass logs
    for _ in 1:5
        step!(integ)
    end
    return integ.u
end

const p = 2.0
const exact = -0.5 * exp(-p * 0.5)  # d/dp exp(-p t) at t = 0.5; Tsit5 at dt = 0.1 is ~2e-7 off

@testset "AD through steps that send progress updates" begin
    logger = Test.TestLogger(min_level = LogLevel(-1))
    with_logger(() -> progress_steps_solution(p), logger)
    # every step passes the gate, so the AD runs below go through the logging path
    @test count(r -> get(r.kwargs, :progress, nothing) isa Float64, logger.logs) == 5

    # (Enzyme forward mode fails on init/step! with or without progress, so it's left out.)
    @testset "$name" for (name, backend) in (
            ("ForwardDiff", AutoForwardDiff()),
            ("Mooncake reverse", AutoMooncake(; config = nothing)),
            ("Mooncake forward", AutoMooncakeForward(; config = nothing)),
            ("Enzyme reverse", AutoEnzyme(; mode = Enzyme.Reverse)),
        )
        logger = Test.TestLogger(min_level = LogLevel(-1))
        d = with_logger(() -> DI.derivative(progress_steps_solution, backend, p), logger)
        @test d ≈ exact rtol = 1.0e-5
        # the intermediate updates report a plain Float64 also when t carries AD types
        # (the start of the solve reports the integer 0, the end "done")
        progress = [r.kwargs[:progress] for r in logger.logs if haskey(r.kwargs, :progress)]
        intermediate = filter(x -> x isa Real && x !== 0, progress)
        @test !isempty(intermediate) && all(x -> x isa Float64, intermediate)
    end
end

# Dual numbers in the time itself (derivative with respect to the end of tspan): the
# fraction is computed on dual times and must still reach the logger as a Float64.
function progress_dual_time(tf)
    prob = ODEProblem((u, p, t) -> -u, 1.0, (zero(tf), tf))
    integ = init(
        prob, Tsit5(); progress = true, save_everystep = false, dt = tf / 10, adaptive = false
    )
    integ.opts.progress_min_interval_ns = 0
    solve!(integ)
    return integ.u
end

@testset "ForwardDiff with dual time" begin
    logger = Test.TestLogger(min_level = LogLevel(-1))
    d = with_logger(() -> ForwardDiff.derivative(progress_dual_time, 1.0), logger)
    @test d ≈ -exp(-1.0) rtol = 1.0e-5
    progress = [r.kwargs[:progress] for r in logger.logs if haskey(r.kwargs, :progress)]
    intermediate = filter(x -> x isa Real && x !== 0, progress)
    @test length(intermediate) >= 9 && all(x -> x isa Float64, intermediate)
end
