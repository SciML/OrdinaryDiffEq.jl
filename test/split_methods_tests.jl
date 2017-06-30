using DiffEqBase, OrdinaryDiffEq, Base.Test

f1 = (t,u) -> 2u
f2 = (t,u) -> 2u

prob = ODEProblem((f1,f2),1.0,(0.0,1.0))

sol = solve(prob,SplitEuler(),dt=1/10)

f3 = (t,u) -> 4u
prob2 = ODEProblem(f3,1.0,(0.0,1.0))
sol2 = solve(prob2,Euler(),dt=1/10)

@test sol2[end] == sol[end]
@test sol2(0.345) == sol(0.345)

u = rand(4,2)
f1 = (t,u,du) -> du.=2u
f2 = (t,u,du) -> du.=2u
prob = ODEProblem((f1,f2),u,(0.0,1.0))
sol = solve(prob,SplitEuler(),dt=1/10)

f3 = (t,u,du) -> du.=4u
prob2 = ODEProblem(f3,u,(0.0,1.0))
sol2 = solve(prob2,Euler(),dt=1/10)

@test sol2[end] == sol[end]
@test sol2(0.345) == sol(0.345)
