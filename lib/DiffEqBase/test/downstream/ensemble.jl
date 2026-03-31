using StochasticDiffEq, DiffEqBase, OrdinaryDiffEq
using Test, Random

import SDEProblemLibrary: prob_sde_2Dlinear, prob_sde_additivesystem, prob_sde_lorenz
import ODEProblemLibrary: prob_ode_linear

prob = prob_sde_2Dlinear
prob2 = EnsembleProblem(prob)
sim = solve(prob2, SRIW1(), dt = 1 // 2^(3), trajectories = 10)

@test sim.u[1] isa DiffEqBase.RODESolution
@test sim[1, 2] isa Matrix
@test sim[1, 2, 1] isa Float64

sim = solve(prob2, SRIW1(), EnsembleThreads(), dt = 1 // 2^(3), trajectories = 10)
err_sim = DiffEqBase.calculate_ensemble_errors(sim; weak_dense_errors = true)
@test length(sim) == 10

sim = solve(
    prob2, SRIW1(), EnsembleThreads(), dt = 1 // 2^(3), trajectories = 10,
    batch_size = 2
)
err_sim = DiffEqBase.calculate_ensemble_errors(sim; weak_dense_errors = true)
@test length(sim) == 10

sim = solve(
    prob2, SRIW1(), EnsembleThreads(), dt = 1 // 2^(3), adaptive = false,
    trajectories = 10
)
err_sim = DiffEqBase.calculate_ensemble_errors(sim; weak_timeseries_errors = true)

sim = solve(prob2, SRIW1(), EnsembleThreads(), dt = 1 // 2^(3), trajectories = 10)
DiffEqBase.calculate_ensemble_errors(sim)
@test length(sim) == 10

sim = solve(prob2, SRIW1(), EnsembleSplitThreads(), dt = 1 // 2^(3), trajectories = 10)
DiffEqBase.calculate_ensemble_errors(sim)
@test length(sim) == 10

sim = solve(prob2, SRIW1(), EnsembleSerial(), dt = 1 // 2^(3), trajectories = 10)
DiffEqBase.calculate_ensemble_errors(sim)
@test length(sim) == 10

prob = prob_sde_additivesystem
prob2 = EnsembleProblem(prob)
sim = solve(prob2, SRA1(), dt = 1 // 2^(3), trajectories = 10)
DiffEqBase.calculate_ensemble_errors(sim)

output_func = function (sol, i)
    return last(last(sol))^2, false
end
prob2 = EnsembleProblem(prob, output_func = output_func)
sim = solve(prob2, SRA1(), dt = 1 // 2^(3), trajectories = 10)

prob = prob_sde_lorenz
prob2 = EnsembleProblem(prob)
sim = solve(prob2, SRIW1(), dt = 1 // 2^(3), trajectories = 10)

output_func = function (sol, i)
    return last(sol), false
end

prob = prob_ode_linear
prob_func = function (prob, i, repeat)
    return ODEProblem(prob.f, rand() * prob.u0, prob.tspan, 1.01)
end

Random.seed!(100)
reduction = function (u, batch, I)
    u = append!(u, batch)
    μ = sum(u) / length(u)
    σ² = sum((x - μ)^2 for x in u) / length(u)
    return u, ((σ² / sqrt(last(I))) / μ < 0.5) ? true : false
end

prob2 = EnsembleProblem(
    prob, prob_func = prob_func, output_func = output_func,
    reduction = reduction, u_init = Vector{Float64}(),
    safetycopy = false
)
sim = solve(prob2, Tsit5(), trajectories = 10000, batch_size = 20)
@test sim.converged == true

prob_func = function (prob, i, repeat)
    return ODEProblem(prob.f, (1 + i / 100) * prob.u0, prob.tspan, 1.01)
end

reduction = function (u, batch, I)
    u = append!(u, batch)
    return u, false
end

prob2 = EnsembleProblem(
    prob, prob_func = prob_func, output_func = output_func,
    reduction = reduction, u_init = Vector{Float64}()
)
sim = solve(prob2, Tsit5(), trajectories = 100, batch_size = 20)
@test sim.converged == false

reduction = function (u, batch, I)
    return u + sum(batch), false
end
prob2 = EnsembleProblem(
    prob, prob_func = prob_func, output_func = output_func,
    reduction = reduction, u_init = 0.0
)
sim2 = solve(prob2, Tsit5(), trajectories = 100, batch_size = 20)
@test sim2.converged == false
@test sum(sim.u) / length(sim.u) ≈ sim2.u / 100

struct SomeUserType end
output_func = function (sol, i)
    return (SomeUserType(), false)
end
prob2 = EnsembleProblem(prob, prob_func = prob_func, output_func = output_func)
sim2 = solve(prob2, Tsit5(), trajectories = 2)
@test sim2.converged && typeof(sim2.u) == Vector{SomeUserType}
