using StochasticDiffEq, Test, Random
using DiffEqGPU
using CUDA

function weak_error(
        prob, alg, numtraj, batchsize, f_true, trange;
        abstol = 1, reltol = 0, ensemblealg = EnsembleCPUArray()
    )
    sol = @time solve(
        prob, alg, ensemblealg,
        dt = 0.01f0, adaptive = true, abstol = abstol, reltol = reltol,
        trajectories = numtraj, batch_size = batchsize,
        saveat = trange
    )
    computed_exp = (sol.u / numtraj)[1, :]
    true_exp = f_true.(trange)
    return sum((computed_exp - true_exp) .^ 2) / length(trange)
end

function prob_func(prob, i, repeat)
    #(i%50000 == 0) && @show i
    return remake(prob, seed = seeds[i])
end

function reduction(u, batch, I)
    return u .+ sum(batch), false
end

function output_func(sol, i)
    #h1(asinh(sol.u[end][1])),false
    return h1.(asinh.(sol)), false
end

# prob 1
u₀ = [0.0f0]
tspan = (0.0f0, 2.0f0)
h1(z) = z^3 - 6 * z^2 + 8 * z

function f1!(du, u, p, t)
    @inbounds begin
        du[1] = 5.0f-1 * u[1] + sqrt(u[1]^2 + 1)
    end
    return nothing
end

function g1!(du, u, p, t)
    @inbounds begin
        du[1] = sqrt(u[1]^2 + 1)
    end
    return nothing
end

f_true1(t) = t^3 - 3 * t^2 + 2 * t

prob1 = SDEProblem(f1!, g1!, u₀, tspan)
ensemble_prob1 = EnsembleProblem(
    prob1;
    output_func = output_func,
    prob_func = prob_func,
    reduction = reduction,
    u_init = Vector{eltype(prob1.u0)}([0.0]),
    safetycopy = false
)

# prob 2
u₀ = [0.1f0, 0.1f0]
function f2!(du, u, p, t)
    @inbounds begin
        du[1] = 1.5f-0 * u[1]
        du[2] = 1.5f-0 * u[2]
    end
    return nothing
end
function g2!(du, u, p, t)
    @inbounds begin
        du[1] = 1.0f-1 * u[1]
        du[2] = 1.0f-1 * u[2]
    end
    return nothing
end

f_true2(t) = 1.0f-1 * exp(1.5f-1 * t) #1//100*exp(301//100*t)

h2(z) = z #z^2

prob2 = SDEProblem(f2!, g2!, u₀, tspan)
ensemble_prob2 = EnsembleProblem(
    prob2;
    output_func = (sol, i) -> (h2.(sol), false),
    prob_func = prob_func,
    reduction = reduction,
    u_init = Vector{eltype(prob2.u0)}([0.0, 0.0]),
    safetycopy = false
)

tsave = 0.0f0:0.05f0:2.0f0

probs = Vector{EnsembleProblem}(undef, 2)
probs[1] = ensemble_prob1
probs[2] = ensemble_prob2

ftrue = Vector{}(undef, 2)
ftrue[1] = f_true1
ftrue[2] = f_true2

numtraj = Int(1.0e6)
seed = 100
Random.seed!(seed)
seeds = rand(UInt, numtraj)

for i in 1:2
    @show i

    err1 = weak_error(
        probs[i], DRI1NM(), numtraj, Int(1.0e4), ftrue[i], tsave, abstol = 1.0f0,
        reltol = 1.0f0, ensemblealg = EnsembleGPUArray(CUDA.CUDABackend())
    )
    @show err1
    # err2 = weak_error(probs[i],DRI1NM(),numtraj,Int(1e4),ftrue[i],tsave,abstol=0.1f0,reltol=0.1f0, ensemblealg=EnsembleGPUArray(CUDA.CUDABackend()))
    # @show err2
    err3 = weak_error(
        probs[i], DRI1NM(), numtraj, Int(1.0e4), ftrue[i], tsave, abstol = 0.01f0,
        reltol = 0.01f0, ensemblealg = EnsembleGPUArray(CUDA.CUDABackend())
    )
    @show err3
    @test err1 > err3
    println("")
end

#
# sol = @time solve(probs[1],DRI1NM(),EnsembleGPUArray(CUDA.CUDABackend()),
#   dt=0.001f0,adaptive=false,abstol=0.1f0,reltol=0.1f0,
#   trajectories=numtraj,batch_size=Int(1e1),
#   saveat = tsave
#   )
# computed_exp = (sol.u/numtraj)[1,:]
# true_exp = f_true1.(tsave)
#
# @show  sum((computed_exp-true_exp).^2)/length(tsave)
#
# using Plots; plot(tsave, true_exp)
# plot!(tsave, computed_exp)
