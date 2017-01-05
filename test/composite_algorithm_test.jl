using OrdinaryDiffEq, DiffEqProblemLibrary
prob = prob_ode_2Dlinear
choice_function(integrator) = (Int(integrator.t<0.5) + 1)
alg_double = CompositeAlgorithm((Tsit5(),Tsit5()),choice_function)
alg_switch = CompositeAlgorithm((Vern7(),Tsit5()),choice_function)

integrator = init(prob,alg_switch,dt=1/8)
for i in take(integrator,4) end
step!(integrator)
step!(integrator)
step!(integrator)

@time sol = solve(prob,alg_double)
@time sol = solve(prob,Tsit5())

using BenchmarkTools

@benchmark sol = solve(prob,alg_double)
@benchmark sol = solve(prob,Tsit5())
