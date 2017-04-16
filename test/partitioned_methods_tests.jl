using DiffEqBase, OrdinaryDiffEq, Base.Test

u0 = 3ones(4)
v0 = ones(4)
f1 = function (t,u,v,du)
  du .= v
end
f2 = function (t,u,v,dv)
  dv .= -2u
end

prob = PartitionedODEProblem((f1,f2),(u0,v0),(0.0,5.0))

sol = solve(prob,SymplecticEuler(),dt=1/100)

interp_time = 0:0.001:5
interps = sol(interp_time)

prob = SecondOrderODEProblem(f2,u0,v0,(0.0,5.0))
sol2 = solve(prob,SymplecticEuler(),dt=1/100)

@test sol[end][1] == sol2[end][1]
@test sol[end][5] == sol2[end][5]
