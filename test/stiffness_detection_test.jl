using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

prob = prob_ode_vanderpol

choice_function(integrator) = (Int(integrator.dt<0.001) + 1)
alg_switch = CompositeAlgorithm((Tsit5(),Rodas5()),choice_function)
sol = solve(prob, alg_switch, reltol=1e-3, abstol=1e-3);
@test length(sol.t) < 10
