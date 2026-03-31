using OrdinaryDiffEq, Test
f = (u, p, t) -> u * p[1]
prob = ODEProblem(f, 1.01, (0.0, 1.0))
@test_throws SciMLBase.NullParameterIndexError solve(prob, Tsit5())
