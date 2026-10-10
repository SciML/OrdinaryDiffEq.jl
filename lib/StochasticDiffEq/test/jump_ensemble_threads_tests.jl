using StochasticDiffEq, JumpProcesses, SciMLBase, Test

# A seeded jump-diffusion ensemble must give the same trajectories serially and with
# threads, under the default jump-state policy: reuse on the main thread, copy on others.
@test Threads.nthreads() > 1

f!(du, u, p, t) = (du .= 0; nothing)
g!(du, u, p, t) = (du .= 0.1; nothing)
birth_rate(u, p, t) = p[1]
birth!(integrator) = (integrator.u[1] += 1; nothing)

const T_END = 2.0
const TRAJECTORIES = 64

jprob = JumpProblem(
    SDEProblem(f!, g!, [0.0], (0.0, T_END), [5.0]), Direct(),
    ConstantRateJump(birth_rate, birth!)
)

function check_serial_matches_threads(ensemble; kwargs...)
    serial = solve(
        ensemble, EM(), EnsembleSerial(); dt = 0.01, trajectories = TRAJECTORIES, kwargs...
    )
    threaded = solve(
        ensemble, EM(), EnsembleThreads(); dt = 0.01, trajectories = TRAJECTORIES,
        kwargs...
    )
    for sols in (serial, threaded)
        @test all(SciMLBase.successful_retcode, sols.u)
        @test all(sol -> sol.t[end] == T_END, sols.u)
    end
    @test all(i -> serial.u[i].t == threaded.u[i].t, 1:TRAJECTORIES)
    @test all(i -> serial.u[i].u == threaded.u[i].u, 1:TRAJECTORIES)
    # Each trajectory follows its own path.
    @test allunique(sol.u[end] for sol in serial.u)
    return nothing
end

@testset "Seeded jump-diffusion ensembles agree serially and with threads" begin
    # Each trajectory gets its own seed, through the wrapped problem. A custom `prob_func`
    # turns on SciMLBase's per-trajectory safety copies.
    seeded(prob, ctx) = remake(prob; prob = remake(prob.prob; seed = 1000 + ctx.sim_id))
    check_serial_matches_threads(EnsembleProblem(jprob; prob_func = seeded))
end

@testset "Default ensembles reuse jump problems and agree serially and with threads" begin
    # With the default `prob_func` there are no per-trajectory safety copies: a serial
    # solve reuses the problem for every trajectory, and a threaded solve reuses one copy
    # per task. The ensemble seed gives each trajectory its own seed.
    ensemble = EnsembleProblem(jprob)
    @test !ensemble.safetycopy
    check_serial_matches_threads(ensemble; seed = 1234)
end
