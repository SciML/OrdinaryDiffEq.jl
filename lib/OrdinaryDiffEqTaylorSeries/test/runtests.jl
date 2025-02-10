using OrdinaryDiffEqTaylorSeries

f(u, p, t) = cos(u)
f!(du, u, p, t) = du .= cos.(u)
odef = ODEFunction{false, SciMLBase.NoSpecialize}(f)
odef! = ODEFunction{true, SciMLBase.NoSpecialize}(f!)

u0 = 0.0
u0! = [0.0]
prob = ODEProblem(odef, u0, (0.0, 10.0))
prob! = ODEProblem(odef!, u0!, (0.0, 10.0))
sol = solve(prob, ExplicitTaylor2(), dt=0.01)
sol! = solve(prob!, ExplicitTaylor2(), dt=0.01)
