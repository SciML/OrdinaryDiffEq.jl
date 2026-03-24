using OrdinaryDiffEq
function f(du, u, p, t)
    return du[1] = 1.01 * u[1]
end
u0 = [0.0, 0.0]
tspan = (0.0, 1.0)
prob = ODEProblem(f, u0, tspan)
n = 100

initial_conditions = range(0, stop = 1, length = n)
function prob_func(prob, i, repeat)
    prob.u0[1] = initial_conditions[i]
    return prob
end
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func)
sim_1 = solve(
    ensemble_prob, Tsit5(), EnsembleThreads(),
    trajectories = 100
)
sim_2 = solve(
    ensemble_prob, Tsit5(), EnsembleDistributed(),
    trajectories = 100
)
ss_sol_1 = hcat(collect(EnsembleAnalysis.get_timepoint(sim_1, 0))...);
ss_sol_2 = hcat(collect(EnsembleAnalysis.get_timepoint(sim_2, 0))...);

ss_sol_1[1, :] == initial_conditions
ss_sol_2[1, :] == initial_conditions
