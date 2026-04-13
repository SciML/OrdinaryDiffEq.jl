using DiffEqGPU, StochasticDiffEq, Test, DiffEqNoiseProcess
using CUDA
using Random

function f!(du, u, p, t)
    @inbounds begin
        du[1] = p[1] * u[1]
        du[2] = -p[1] * u[2]
    end
    return nothing
end

function scalar_noise!(du, u, p, t)
    @inbounds begin
        du[1] = p[2]
        du[2] = -p[2]
    end
    return nothing
end

function prob_func(prob, i, repeat)
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

CUDA.allowscalar(false)
u0 = [-0.1f0, 0.0f0]
tspan = (0.0f0, 100.0f0)
p = [1.9f0, 0.1f0]

prob = SDEProblem(f!, scalar_noise!, u0, tspan, p, noise = W)
ensembleprob = EnsembleProblem(prob, prob_func = prob_func)

@info "scalar SDE"

#Performance check with nvvp
# CUDAnative.CUDAdrv.@profile
# check either on CPU with EnsembleCPUArray() or on GPU with EnsembleGPUArray()
#@test_nowarn sol = @time solve(ensembleprob,DRI1(),EnsembleCPUArray(),trajectories=numtraj)
sol = @time solve(ensembleprob, DRI1NM(), EnsembleGPUArray(CUDA.CUDABackend()), trajectories = numtraj)
