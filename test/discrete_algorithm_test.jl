using OrdinaryDiffEq, DiffEqProblemLibrary, Base.Test

sol =solve(prob_ode_linear,Discrete())

@test sol[1] == sol[end]

@test sol(0.5:0.1:0.7) == [sol[1],sol[1],sol[1]]

sol =solve(prob_ode_2Dlinear,Discrete())

@test sol[1] == sol[end]

@test sol(0.5) == sol[1]

sol =solve(prob_ode_linear,FunctionMap())

@test sol[end] ≈ .505

sol =solve(prob_ode_linear,FunctionMap(scale_by_time=true),dt=1/4)

sol2 =solve(prob_ode_linear,Euler(),dt=1/4)

@test sol[end] == sol2[end]

@test sol(0.53) != sol2(0.53)

sol =solve(prob_ode_2Dlinear,FunctionMap())

sol2 =solve(prob_ode_2Dlinear,Euler(),dt=1)

sol =solve(prob_ode_2Dlinear,FunctionMap(scale_by_time=true),dt=1/4)

sol2 =solve(prob_ode_2Dlinear,Euler(),dt=1/4)

@test sol[end] == sol2[end]

@test sol(0.35) != sol2(0.53)
