using Distributed
addprocs(2)
println("There are $(nprocs()) processes")

@everywhere begin
    using Pkg
    Pkg.activate("downstream")
    Pkg.develop(PackageSpec(path = joinpath(pwd(), "..")))
    Pkg.instantiate()
    using OrdinaryDiffEq
    prob = ODEProblem((u, p, t) -> 1.01u, 0.5, (0.0, 1.0))
    u0s = [rand() * prob.u0 for i in 1:2]
    function simple_prob_func(prob, i, repeat)
        println("Running trajectory $i")
        ODEProblem(prob.f, u0s[i], prob.tspan)
    end
end

ensemble_prob = EnsembleProblem(prob, prob_func = simple_prob_func)
sim = solve(ensemble_prob, Tsit5(), EnsembleSplitThreads(), trajectories = 2)

@everywhere function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0, 0.0, 0.0]
tspan = (0.0, 100.0)
p = [1, 2.0, 3]
prob = ODEProblem(lorenz!, u0, tspan, p)

@everywhere function lorenz_prob_func(prob, i, repeat)
    prob = remake(prob, tspan = (rand(), 100.0), p = rand(3))
    return prob
end

ensemble_prob = EnsembleProblem(prob, prob_func = lorenz_prob_func, safetycopy = true)

println("Running EnsembleSerial()")
@test length(solve(ensemble_prob, Tsit5(), EnsembleSerial(), trajectories = 100)) == 100
println("Running EnsembleThreads()")
@test length(solve(ensemble_prob, Tsit5(), EnsembleThreads(), trajectories = 100)) == 100
println("Running EnsembleDistributed()")
@test length(solve(ensemble_prob, Tsit5(), EnsembleDistributed(), trajectories = 100)) ==
    100
println("Running EnsembleSplitThreads()")
@test length(solve(ensemble_prob, Tsit5(), EnsembleSplitThreads(), trajectories = 100)) ==
    100
