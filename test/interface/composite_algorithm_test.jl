using OrdinaryDiffEq, Test
using DiffEqProblemLibrary.ODEProblemLibrary: importodeproblems; importodeproblems()
import DiffEqProblemLibrary.ODEProblemLibrary: prob_ode_linear, prob_ode_2Dlinear
prob = prob_ode_2Dlinear
choice_function(integrator) = (Int(integrator.t<0.5) + 1)
alg_double = CompositeAlgorithm((Tsit5(),Tsit5()),choice_function)
alg_double2 = CompositeAlgorithm((Vern6(),Vern6()),choice_function)
alg_switch = CompositeAlgorithm((Tsit5(),Vern7()),choice_function)

@time sol1 = solve(prob_ode_linear,alg_double)
@time sol2 = solve(prob_ode_linear,Tsit5())
@test sol1.t == sol2.t
@test sol1(0.8) == sol2(0.8)

integrator1 = init(prob,alg_double2)
integrator2 = init(prob,Vern6())
solve!(integrator1)
solve!(integrator2)

@test integrator1.sol.t == integrator2.sol.t

sol = solve(prob,alg_switch)
@inferred OrdinaryDiffEq.DiffEqBase.__solve(prob,alg_switch)
