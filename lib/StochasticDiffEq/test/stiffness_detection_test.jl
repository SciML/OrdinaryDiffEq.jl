using StochasticDiffEq, Test, Random
using SDEProblemLibrary: prob_sde_stiffquadito

Random.seed!(100)
prob = prob_sde_stiffquadito
prob = remake(prob; p = (1.0e5, 2.0))
alg = AutoSOSRA2(SKenCarp(), maxstiffstep = 5, maxnonstiffstep = 2, stiffalgfirst = false)
@test StochasticDiffEq.isadaptive(prob, alg)
@test StochasticDiffEq.isadaptive(alg)
@time sol = solve(prob, alg)
@test alg.algs[sol.alg_choice[end]] isa SKenCarp
@test length(unique(sol.alg_choice)) == 2

Random.seed!(100)
prob = prob_sde_stiffquadito
prob = remake(prob; p = (1.0e5, 2.0))
@time sol = solve(
    prob, AutoSOSRI2(ImplicitRKMil(), maxstiffstep = 1, maxnonstiffstep = 10, stiffalgfirst = false)
)
@test length(unique(sol.alg_choice)) == 2
