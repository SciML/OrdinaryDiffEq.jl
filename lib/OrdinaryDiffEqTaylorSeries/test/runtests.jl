using OrdinaryDiffEqTaylorSeries

f(u, p, t) = cos(u)
f!(du, u, p, t) = du .= cos.(u)

u0 = 0.0
u0! = [0.0]
prob = ODEProblem{false, SciMLBase.NoSpecialize}(f, u0, (0.0, 10.0))
prob! = ODEProblem{true, SciMLBase.NoSpecialize}(f!, u0!, (0.0, 10.0))
sol = solve(prob, ExplicitTaylor2(), dt=0.01)
sol! = solve(prob!, ExplicitTaylor2(), dt=0.01)
sol = solve(prob, ExplicitTaylor(order=Val(2)), dt=0.01)
sol! = solve(prob!, ExplicitTaylor(order=Val(2)), dt=0.01)
