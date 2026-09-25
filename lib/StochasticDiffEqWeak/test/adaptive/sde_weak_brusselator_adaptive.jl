using StochasticDiffEqWeak, Test, DiffEqNoiseProcess
using SciMLBase
using Random

function brusselator_f!(du, u, p, t)
    @inbounds begin
        du[1] = (p[1] - 1) * u[1] + p[1] * u[1]^2 + (u[1] + 1)^2 * u[2]
        du[2] = -p[1] * u[1] - p[1] * u[1]^2 - (u[1] + 1)^2 * u[2]
    end
    return nothing
end

function scalar_noise!(du, u, p, t)
    @inbounds begin
        du[1] = p[2] * u[1] * (1 + u[1])
        du[2] = -p[2] * u[1] * (1 + u[1])
    end
    return nothing
end

function additive_scalar_noise!(du, u, p, t)
    @inbounds begin
        du[1] = p[2]
        du[2] = -p[2]
    end
    return nothing
end

function prob_func(prob, ctx)
    i = ctx.sim_id; repeat = ctx.repeat
    Random.seed!(seeds[i])
    W = WienerProcess(0.0, 0.0, 0.0)
    return remake(prob, noise = W)
end

# fix seeds
seed = 100
Random.seed!(seed)
numtraj = 100
seeds = rand(UInt, numtraj)
W = WienerProcess(0.0, 0.0, 0.0)

#CUDA.allowscalar(false)
u0 = [-0.1f0, 0.0f0]
tspan = (0.0f0, 100.0f0)
p = [1.9f0, 0.1f0]

prob = SDEProblem(brusselator_f!, scalar_noise!, u0, tspan, p, noise = W)
ensembleprob = EnsembleProblem(prob; prob_func)

@info "Brusselator"

short_prob = remake(prob; tspan = (0.0f0, 1.0f0))
short_additive_prob = SDEProblem(
    brusselator_f!, additive_scalar_noise!, u0, (0.0f0, 1.0f0), p,
    noise = WienerProcess(0.0, 0.0, 0.0)
)

@testset "$(nameof(typeof(alg))) with scalar Wiener process" for (alg, short_test_prob) in
    ((PL1WM(), short_prob), (PL1WMA(), short_additive_prob))
    sol = solve(short_test_prob, alg; dt = 0.1f0, adaptive = false)
    @test SciMLBase.successful_retcode(sol)
end

#Performance check with nvvp
# CUDAnative.CUDAdrv.@profile
# check either on CPU with EnsembleCPUArray() or on GPU with EnsembleGPUArray()
@test_nowarn sol = @time solve(ensembleprob, DRI1(), EnsembleThreads(), trajectories = numtraj)
#sol = @time solve(ensembleprob,DRI1(),EnsembleGPUArray(),trajectories=numtraj)
