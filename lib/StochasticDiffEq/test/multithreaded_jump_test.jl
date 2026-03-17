using JumpProcesses, StochasticDiffEq, Test

# Requires threads to be effective
# https://github.com/SciML/DifferentialEquations.jl/issues/854
using Base.Threads
@test Threads.nthreads() > 1

function testdrift!(du, u, p, t)
    return du[1] = u[1]
end

function testdiffusion!(du, u, p, t)
    return du[1] = u[1]
end

function testrate(u, p, t)
    return 0.32
end

function testaffect!(integrator)
    return integrator.u[1] += integrator.u[1] * randn()
end

testjump = ConstantRateJump(testrate, testaffect!)

test_sdeprob = SDEProblem(
    testdrift!,
    testdiffusion!,
    [1.0],
    (0.0, 1.0)
)

test_jumpprob = JumpProblem(
    test_sdeprob,
    Direct(),
    testjump
)

test_ensemprob = EnsembleProblem(test_jumpprob)

test_ensemsim = solve(
    test_ensemprob,
    EM(),
    dt = 0.01,
    EnsembleThreads();
    trajectories = 1_000_000,
    alias_jump = false
)
